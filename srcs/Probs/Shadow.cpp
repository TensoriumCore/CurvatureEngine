#include "MetricAVX.h"
#include <Geodesics.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <omp.h>
#include <vector>

extern float a;
extern float (*geodesic_points)[5];

using namespace curvatureengine::simd;

struct StageDerivativesVec {
  VEC_TYPE dx[NDIM];
  VEC_TYPE dv[NDIM];
};

struct PixelResult {
  double r[4];
  double g[4];
  double b[4];
};

struct Tetrad {
  double e[4][4];
};

static inline double dot_g(const double g[4][4], const double u[4],
                           const double v[4]) {
  double s = 0.0;
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu)
      s += g[mu][nu] * u[mu] * v[nu];
  return s;
}

static inline void axpy(double v[4], double c, const double u[4]) {
  for (int i = 0; i < 4; ++i)
    v[i] += c * u[i];
}

static inline void axpy_neg(double v[4], double c, const double u[4]) {
  for (int i = 0; i < 4; ++i)
    v[i] -= c * u[i];
}

static inline void project_orthogonal_to_u(const double g[4][4],
                                           const double u[4], double v[4]) {
  const double uv = dot_g(g, u, v);
  axpy(v, uv, u);
}

static inline void normalize_spacelike(const double g[4][4], double v[4]) {
  double vv = dot_g(g, v, v);
  double inv = 1.0 / std::sqrt(std::max(vv, 1e-300));
  for (int i = 0; i < 4; ++i)
    v[i] *= inv;
}

static inline double clamp255(double x) {
  if (x < 0.0)
    return 0.0;
  if (x > 255.0)
    return 255.0;
  return x;
}

static inline void build_camera_tetrad_zamo(const double g[4][4],
                                            const double ginv[4][4],
                                            Tetrad &T) {
  const double gtt_inv = ginv[0][0];
  const double alpha = 1.0 / std::sqrt(std::max(-gtt_inv, 1e-300));

  const double beta_r = ginv[0][1] / gtt_inv;
  const double beta_th = ginv[0][2] / gtt_inv;
  const double beta_ph = ginv[0][3] / gtt_inv;

  double u[4] = {1.0 / alpha, -beta_r / alpha, -beta_th / alpha,
                 -beta_ph / alpha};

  double uu = dot_g(g, u, u);
  if (uu < 0.0) {
    double s = 1.0 / std::sqrt(-uu);
    for (int i = 0; i < 4; ++i)
      u[i] *= s;
  }

  for (int mu = 0; mu < 4; ++mu)
    T.e[0][mu] = u[mu];

  double er[4] = {0.0, 1.0, 0.0, 0.0};
  double eth[4] = {0.0, 0.0, 1.0, 0.0};
  double eph[4] = {0.0, 0.0, 0.0, 1.0};

  project_orthogonal_to_u(g, u, er);
  normalize_spacelike(g, er);
  project_orthogonal_to_u(g, u, eth);
  {
    double c = dot_g(g, er, eth);
    axpy_neg(eth, c, er);
    normalize_spacelike(g, eth);
  }

  project_orthogonal_to_u(g, u, eph);
  {
    double c1 = dot_g(g, er, eph);
    double c2 = dot_g(g, eth, eph);
    axpy_neg(eph, c1, er);
    axpy_neg(eph, c2, eth);
    normalize_spacelike(g, eph);
  }

  for (int mu = 0; mu < 4; ++mu) {
    T.e[1][mu] = er[mu];
    T.e[2][mu] = eth[mu];
    T.e[3][mu] = eph[mu];
  }
}

static inline void make_camera_rays_4(const Tetrad &T, const double sx[4],
                                      const double sy[4], double fov,
                                      double k_out[4][4]) {
  const double tan_half_fov = std::tan(0.5 * fov);

  for (int i = 0; i < 4; ++i) {
    double dx = sx[i] * tan_half_fov; 
    double dy = sy[i] * tan_half_fov;

    double norm = std::sqrt(1.0 + dx * dx + dy * dy);
    double inv_norm = 1.0 / norm;

    double n_fwd = 1.0 * inv_norm;
    double n_right = dx * inv_norm;
    double n_up = dy * inv_norm;

    for (int mu = 0; mu < 4; ++mu) {
      k_out[i][mu] = T.e[0][mu] - n_fwd * T.e[1][mu] + n_right * T.e[3][mu] -
                     n_up * T.e[2][mu];
    }
  }
}

static void compute_vector_stage_local(const VEC_TYPE state_x[NDIM],
                                       const VEC_TYPE state_v[NDIM],
                                       ChristoffelEvalFnVec eval_fn, void *ctx,
                                       StageDerivativesVec &stage) {
  ChristoffelTensorVec gamma;
  if (eval_fn)
    eval_fn(state_x, gamma, ctx);

  for (int mu = 0; mu < NDIM; ++mu) {
    stage.dx[mu] = state_v[mu];
    VEC_TYPE acc = zero();
    if (eval_fn) {
      for (int alpha = 0; alpha < NDIM; ++alpha) {
        VEC_TYPE v_alpha = state_v[alpha];
        for (int beta = 0; beta < NDIM; ++beta) {
          acc = fnmadd(gamma[mu][alpha][beta] * v_alpha, state_v[beta], acc);
        }
      }
    }
    stage.dv[mu] = acc;
  }
}

static inline double compute_risco(double a) {
  double aa = std::min(0.999999999, std::max(-0.999999999, a));
  double Z1 = 1.0 + std::cbrt(1.0 - aa * aa) *
                        (std::cbrt(1.0 + aa) + std::cbrt(1.0 - aa));
  double Z2 = std::sqrt(3.0 * aa * aa + Z1 * Z1);
  double s = (aa >= 0.0) ? 1.0 : -1.0;
  return 3.0 + Z2 - s * std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2));
}

static inline double nt_flux_simple(double r, double risco) {
  double inv_r = 1.0 / std::max(r, 1e-12);
  double f = 1.0 - std::sqrt(risco * inv_r);
  if (f <= 0.0)
    return 0.0;
  return f * inv_r * inv_r * inv_r;
}

static inline double safe_pos(double x) { return (x > 1e-30) ? x : 1e-30; }

static inline double clamp01(double x) {
  if (x < 0.0)
    return 0.0;
  if (x > 1.0)
    return 1.0;
  return x;
}

static inline void tone_map_rgb(double Ir, double Ig, double Ib, double &out_r,
                                double &out_g, double &out_b) {
  auto TM = [](double x) {
    x = std::max(0.0, x);
    double y = x / (1.0 + x);
    return 255.0 * y;
  };
  out_r = clamp255(TM(Ir));
  out_g = clamp255(TM(Ig));
  out_b = clamp255(TM(Ib));
}

static inline void blackbody_rgb(double T, double &Red, double &Green,
                                 double &Blue) {
  double t = T / 100.0;
  double r, g, b;

  if (t <= 66.0)
    r = 1.0;
  else {
    double x = t - 60.0;
    r = 1.292936186062745 * std::pow(x, -0.1332047592);
  }

  if (t <= 66.0) {
    double x = std::max(t, 1.0);
    g = 0.3900815787690196 * std::log(x) - 0.6318414437886275;
  } else {
    double x = t - 60.0;
    g = 1.129890860895294 * std::pow(x, -0.0755148492);
  }

  if (t >= 66.0)
    b = 1.0;
  else if (t <= 19.0)
    b = 0.0;
  else {
    double x = t - 10.0;
    b = 0.5432067891101961 * std::log(x) - 1.19625408914;
  }

  Red = clamp01(r);
  Green = clamp01(g);
  Blue = clamp01(b);
}

static inline double model_ne(double r, double theta) {
  const double PI_HALF = 1.5707963267948966;
  const double n0 = 2.0;
  const double p = 1.15;
  const double sigma0 = 0.045;
  double dth = theta - PI_HALF;
  double sigma = sigma0 * (1.0 + 0.6 * std::exp(-(r - 6.0) / 18.0));
  double wth = std::exp(-0.5 * (dth / sigma) * (dth / sigma));
  return n0 * std::pow(std::max(r, 1.0), -p) * wth;
}

static inline double model_Te(double r) {
  const double T0 = 2.2e4;
  const double q = 0.60;
  return T0 * std::pow(std::max(r, 1.0), -q);
}

static inline double model_B(double r, double theta, double ne) {
  const double PI_HALF = 1.5707963267948966;
  const double B0 = 0.18;
  const double m = 1.05;
  double dth = std::abs(theta - PI_HALF);
  double w = std::exp(-dth / 0.20);
  double br = B0 * std::pow(std::max(r, 1.0), -m) * w;
  double bn = std::pow(std::max(ne, 1e-12), 0.5);
  return br * bn;
}

static inline double model_Thetae(double Te) {
  const double inv = 1.0 / 6.0e5;
  return std::max(1e-4, Te * inv);
}

static inline double sync_kernel(double x) {
  if (x < 1e-8)
    return 0.0;
  double a = std::pow(x, 1.0 / 3.0);
  double b = std::exp(-x);
  return a * b;
}

static inline void emiss_abs_synch_approx(double ne, double Thetae, double B,
                                          double nu, double &jnu, double &anu) {
  double nu_c = safe_pos(30.0 * Thetae * Thetae * B);
  double x = nu / nu_c;
  double K = sync_kernel(x);

  const double j0 = 6.0e-7;
  jnu = j0 * ne * B * K;

  const double a0 = 8.0e-3;
  anu = a0 * ne * B / safe_pos(nu_c) * K;
}

static inline double compute_Kepler_Omega(double r, double a_d) {
  return 1.0 / (std::pow(r, 1.5) + a_d);
}

static inline bool fluid_ut_uph_from_metric(double g_tt, double g_tph,
                                            double g_phph, double Omega,
                                            double &ut, double &uph,
                                            double &u_cov_t, double &u_cov_ph) {
  double norm_u = -(g_tt + 2.0 * g_tph * Omega + g_phph * Omega * Omega);
  if (!(norm_u > 0.0))
    return false;
  ut = 1.0 / std::sqrt(norm_u);
  uph = Omega * ut;
  u_cov_t = g_tt * ut + g_tph * uph;
  u_cov_ph = g_tph * ut + g_phph * uph;
  return true;
}

static inline void rt_step_exact(double dl, double Jinv, double Ainv,
                                 double &Ibar, double &tau) {
  double dtau = Ainv * dl;
  double att = std::exp(-dtau);
  if (Ainv > 1e-14) {
    double S = Jinv / Ainv;
    Ibar = Ibar * att + S * (1.0 - att);
  } else {
    Ibar = Ibar + Jinv * dl;
  }
  tau += dtau;
}

static inline void accumulate_grrt_step(const VEC_TYPE x_prev[NDIM],
                                        const VEC_TYPE x_curr[NDIM],
                                        const VEC_TYPE k_curr[NDIM],
                                        const double nu_obs[4], double dlambda,
                                        double Ibar[4], double tau[4],
                                        double WTe[4], double Wsum[4]) {
  const double a_d = (double)a;
  double risco = compute_risco(a_d);

  alignas(32) double rP[4], rC[4], thP[4], thC[4];
  store(rP, x_prev[1]);
  store(rC, x_curr[1]);
  store(thP, x_prev[2]);
  store(thC, x_curr[2]);

  alignas(32) double kt[4], kr[4], kth[4], kph[4];
  store(kt, k_curr[0]);
  store(kr, k_curr[1]);
  store(kth, k_curr[2]);
  store(kph, k_curr[3]);

  alignas(32) double rM[4], thM[4];
  for (int i = 0; i < 4; ++i) {
    rM[i] = 0.5 * (rP[i] + rC[i]);
    thM[i] = 0.5 * (thP[i] + thC[i]);
  }

  alignas(32) double Xt[4] = {0}, Xph[4] = {0};
  alignas(32) double Xr[4] = {rM[0], rM[1], rM[2], rM[3]};
  alignas(32) double Xth[4] = {thM[0], thM[1], thM[2], thM[3]};

  VEC_TYPE X_mid[NDIM] = {load(Xt), load(Xr), load(Xth), load(Xph)};
  VEC_TYPE gS[NDIM][NDIM], ginvS[NDIM][NDIM];
  MetricSIMD::calculate_metric_avx(X_mid, gS, ginvS);

  alignas(32) double gtt[4], gtph[4], gphph[4];
  alignas(32) double gtr[4], grr[4], grph[4];
  alignas(32) double gtht[4], gthr[4], gthph[4];

  store(gtt, gS[0][0]);
  store(gtph, gS[0][3]);
  store(gphph, gS[3][3]);
  store(gtr, gS[0][1]);
  store(grr, gS[1][1]);
  store(grph, gS[1][3]);
  store(gtht, gS[2][0]);
  store(gthr, gS[2][1]);
  store(gthph, gS[2][3]);

  const double sigma_th = 0.03;
  const double T_max = 225000.0;
  const double inflow_strength = 0.25;

  const double kappa0 = 1.05;
  const double emiss_scale = 2.5e-12;

  for (int i = 0; i < 4; ++i) {
    double r = rM[i];
    if (r < risco || r > 40.0)
      continue;

    double dth = thM[i] - 1.57079632679;
    double density_profile =
        std::exp(-0.5 * (dth / sigma_th) * (dth / sigma_th));
    if (density_profile < 1e-5)
      continue;

    double flux_nt = nt_flux_simple(r, risco);
    if (flux_nt <= 0.0)
      continue;

    double Te = T_max * std::pow(flux_nt, 0.25);
    double Omega = 1.0 / (std::pow(r, 1.5) + a_d);

    double f = std::max(0.0, 1.0 - risco / r);
    double ur0 = -inflow_strength * f * std::sqrt(2.0 / r);

    double A1 = gtt[i] + 2.0 * gtph[i] * Omega + gphph[i] * Omega * Omega;
    double B1 = 2.0 * ur0 * (gtr[i] + grph[i] * Omega);
    double C1 = grr[i] * ur0 * ur0;
    double norm = -(A1 + B1 + C1);
    if (norm <= 0.0)
      continue;

    double ut = 1.0 / std::sqrt(norm);
    double uph = Omega * ut;
    double ur = ur0 * ut;

    double uct = gtt[i] * ut + gtph[i] * uph + gtr[i] * ur;
    double ucr = gtr[i] * ut + grr[i] * ur + grph[i] * uph;
    double ucth = gtht[i] * ut + gthr[i] * ur + gthph[i] * uph;
    double ucph = gtph[i] * ut + gphph[i] * uph + grph[i] * ur;

    double nu_em = -(uct * kt[i] + ucr * kr[i] + ucth * kth[i] + ucph * kph[i]);
    if (nu_em <= 1e-12)
      continue;

    double redshift = nu_obs[i] / nu_em;
    double T_obs = Te * redshift;

    double anu = kappa0 * density_profile;
    double Snu = emiss_scale * std::pow(Te, 3.0);
    double jnu = anu * Snu;

    double Jinv = jnu / (nu_em * nu_em);
    double Ainv = anu * nu_em;

    double I0 = Ibar[i];
    double tau0 = tau[i];

    rt_step_exact(dlambda, Jinv, Ainv, Ibar[i], tau[i]);

    double dI = Ibar[i] - I0;
    if (dI > 0.0) {
      double w = dI * std::exp(-tau0);
      WTe[i] += T_obs * w;
      Wsum[i] += w;
    }
  }
}

static PixelResult
geodesic_raytrace_physical(VEC_TYPE x[NDIM], VEC_TYPE v[NDIM],
                           const double uk_obs[4], float lambda_max,
                           ChristoffelEvalFnVec eval_fn, void *ctx,
                           double r_horizon, double r_escape) {
  double lambda_lane[4] = {0.0, 0.0, 0.0, 0.0};
  const double l_max = (double)lambda_max;

  double nu_obs[4];
  for (int i = 0; i < 4; ++i)
    nu_obs[i] = safe_pos(-uk_obs[i]);

  double Ibar[4] = {0.0, 0.0, 0.0, 0.0};
  double tau[4] = {0.0, 0.0, 0.0, 0.0};

  double WTe[4] = {0.0, 0.0, 0.0, 0.0};
  double Wsum[4] = {0.0, 0.0, 0.0, 0.0};

  alignas(32) double r_now[4];
  bool finished[4] = {false, false, false, false};
  int alive = 4;

  VEC_TYPE prev_x[NDIM];
  for (int i = 0; i < NDIM; ++i)
    prev_x[i] = x[i];

  while (alive > 0) {
    store(r_now, x[1]);

    for (int i = 0; i < 4; ++i) {
      if (finished[i])
        continue;
      if (std::isnan(r_now[i]) || r_now[i] < r_horizon) {
        finished[i] = true;
        alive--;
        continue;
      }
      if (r_now[i] > r_escape) {
        finished[i] = true;
        alive--;
        continue;
      }
      if (lambda_lane[i] >= l_max) {
        finished[i] = true;
        alive--;
        continue;
      }
      if (tau[i] > 15.0) {
        finished[i] = true;
        alive--;
        continue;
      }
    }

    if (alive == 0)
      break;

    double min_r = 1e100;
    for (int i = 0; i < 4; ++i)
      if (!finished[i])
        min_r = std::min(min_r, r_now[i]);

    double step_val = 0.6;
    if (min_r < 25.0)
      step_val = 0.18;
    if (min_r < 10.0)
      step_val = 0.06;
    if (min_r < 5.0)
      step_val = 0.02;
    if (min_r < 3.0)
      step_val = 0.01;
    if (min_r < r_horizon * 1.5)
      step_val = 0.005;

    double min_rem = 1e300;
    for (int i = 0; i < 4; ++i)
      if (!finished[i])
        min_rem = std::min(min_rem, l_max - lambda_lane[i]);
    step_val = std::min(step_val, min_rem);
    if (!(step_val > 0.0))
      break;

    VEC_TYPE h = broadcast(step_val);

    StageDerivativesVec k1, k2, k3, k4;
    VEC_TYPE tmp_x[NDIM], tmp_v[NDIM];

    compute_vector_stage_local(x, v, eval_fn, ctx, k1);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k1.dx[i] * h * 0.5;
      tmp_v[i] = v[i] + k1.dv[i] * h * 0.5;
    }

    compute_vector_stage_local(tmp_x, tmp_v, eval_fn, ctx, k2);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k2.dx[i] * h * 0.5;
      tmp_v[i] = v[i] + k2.dv[i] * h * 0.5;
    }

    compute_vector_stage_local(tmp_x, tmp_v, eval_fn, ctx, k3);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k3.dx[i] * h;
      tmp_v[i] = v[i] + k3.dv[i] * h;
    }

    compute_vector_stage_local(tmp_x, tmp_v, eval_fn, ctx, k4);

    VEC_TYPE h6 = h * (1.0 / 6.0);

    for (int i = 0; i < NDIM; ++i)
      prev_x[i] = x[i];

    for (int i = 0; i < NDIM; ++i) {
      x[i] =
          x[i] + h6 * (k1.dx[i] + 2.0 * k2.dx[i] + 2.0 * k3.dx[i] + k4.dx[i]);
      v[i] =
          v[i] + h6 * (k1.dv[i] + 2.0 * k2.dv[i] + 2.0 * k3.dv[i] + k4.dv[i]);
    }

    for (int i = 0; i < 4; ++i)
      if (!finished[i])
        lambda_lane[i] += step_val;

    accumulate_grrt_step(prev_x, x, v, nu_obs, step_val, Ibar, tau, WTe, Wsum);
  }

  PixelResult res;
  for (int i = 0; i < 4; ++i) {
    double Te_eff = (Wsum[i] > 0.0) ? (WTe[i] / Wsum[i]) : 9000.0;

    double cr, cg, cb;
    blackbody_rgb(std::max(100.0, Te_eff), cr, cg, cb);

    double Iobs = Ibar[i] * nu_obs[i] * nu_obs[i] * nu_obs[i];

    // Iobs *= 2.0e6;
    double Ir = Iobs * cr;
    double Ig = Iobs * cg;
    double Ib = Iobs * cb;

    tone_map_rgb(Ir, Ig, Ib, res.r[i], res.g[i], res.b[i]);
  }
  return res;
}

void write_ppm_image(const char *filename, int width, int height,
                     const std::vector<uint8_t> &rgb_buffer) {
  FILE *fp = fopen(filename, "wb");
  if (!fp)
    return;
  fprintf(fp, "P6\n%d %d\n255\n", width, height);
  fwrite(rgb_buffer.data(), 1, rgb_buffer.size(), fp);
  fclose(fp);
}

int shadow_prob() {
  a = 0.935f;

  const int WIDTH = 1000;
  const int HEIGHT = 800;
  std::vector<uint8_t> image_buffer(WIDTH * HEIGHT * 3);

  const double r_cam = 75.0;
  const double theta_cam = 1.2865;
  const double phi_cam = 0.0;

  const double fov = 1.3;
  const double aspect = (double)WIDTH / (double)HEIGHT;

  const double r_horizon = 1.0 + std::sqrt(1.0 - (double)a * (double)a);
  const double r_stop = r_horizon + 0.05;

  const int packets_x = WIDTH / 4;

#pragma omp parallel for schedule(dynamic)
  for (int y = 0; y < HEIGHT; ++y) {
    double sy = 2.0 * ((y + 0.5) / (double)HEIGHT) - 1.0;

    sy = -sy;

    for (int p = 0; p < packets_x; ++p) {
      int x_base = p * 4;

      alignas(32) double sx4[4];
      for (int l = 0; l < 4; ++l) {
        double sx = 2.0 * ((x_base + l + 0.5) / (double)WIDTH) - 1.0;
        sx *= aspect;
        sx4[l] = sx;
      }

      VEC_TYPE X_cam[NDIM] = {broadcast(0.0), broadcast(r_cam),
                              broadcast(theta_cam), broadcast(phi_cam)};
      VEC_TYPE gS[NDIM][NDIM], ginvS[NDIM][NDIM];
      MetricSIMD::calculate_metric_avx(X_cam, gS, ginvS);

      alignas(32) double gsc[4], gisc[4];
      double g0[4][4], gi0[4][4];
      for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
          store(gsc, gS[mu][nu]);
          store(gisc, ginvS[mu][nu]);
          g0[mu][nu] = gsc[0];
          gi0[mu][nu] = gisc[0];
        }
      }

      Tetrad T;
      build_camera_tetrad_zamo(g0, gi0, T);

      double sy4[4] = {sy, sy, sy, sy};
      double k_out[4][4];
      make_camera_rays_4(T, sx4, sy4, fov, k_out);

      alignas(32) double kt[4], kr[4], kth[4], kph[4];
      for (int i = 0; i < 4; ++i) {
        kt[i] = k_out[i][0];
        kr[i] = k_out[i][1];
        kth[i] = k_out[i][2];
        kph[i] = k_out[i][3];
      }
      VEC_TYPE v_avx[NDIM] = {load(kt), load(kr), load(kth), load(kph)};

      double u_obs[4] = {T.e[0][0], T.e[0][1], T.e[0][2], T.e[0][3]};
      double u_obs_cov[4] = {0};
      for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int nu = 0; nu < 4; ++nu)
          s += g0[mu][nu] * u_obs[nu];
        u_obs_cov[mu] = s;
      }

      double uk_obs[4];
      for (int i = 0; i < 4; ++i) {
        uk_obs[i] = u_obs_cov[0] * kt[i] + u_obs_cov[1] * kr[i] +
                    u_obs_cov[2] * kth[i] + u_obs_cov[3] * kph[i];
      }

      static int once = 0;
      if (!once) {
        once = 1;
        for (int i = 0; i < 4; ++i) {
          printf("uk_obs[%d]=%.6e  nu_obs=%.6e\n", i, uk_obs[i],
                 safe_pos(-uk_obs[i]));
        }
      }
      VEC_TYPE X_avx[NDIM] = {broadcast(0.0), broadcast(r_cam),
                              broadcast(theta_cam), broadcast(phi_cam)};

      PixelResult res = geodesic_raytrace_physical(
          X_avx, v_avx, uk_obs, 4000.0, evaluate_christoffel_native_avx,
          nullptr, r_stop, 80.0);

      for (int l = 0; l < 4; ++l) {
        int px_idx = (y * WIDTH + (x_base + l)) * 3;
        image_buffer[px_idx + 0] = (uint8_t)res.r[l];
        image_buffer[px_idx + 1] = (uint8_t)res.g[l];
        image_buffer[px_idx + 2] = (uint8_t)res.b[l];
      }
    }
  }

  write_ppm_image("Output/accretion_disk.ppm", WIDTH, HEIGHT, image_buffer);

  alignas(32) static float geodesic_points_fallback[2048][5];
  if (!geodesic_points)
    geodesic_points = geodesic_points_fallback;

  return 0;
}
