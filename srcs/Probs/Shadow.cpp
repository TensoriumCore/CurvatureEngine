#include "MetricAVX.h"
#include "app/Problems.h"
#include "app/RuntimeState.h"
#include "core/Constants.h"
#include "core/GeodesicIntegrator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <omp.h>
#include <string>
#include <vector>

using namespace curvatureengine::simd;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct StageDerivativesVec {
  VEC_TYPE dx[NDIM];
  VEC_TYPE dv[NDIM];
};

struct PixelResult {
  double linear_r[4];
  double linear_g[4];
  double linear_b[4];
  double intensity[4];
  double optical_depth[4];
  double effective_temperature[4];
};

struct Tetrad {
  double e[4][4];
};

struct ShadowRenderConfig {
  int width = 1920;
  int height = 1080;
  int samples_per_pixel = 4;
  int frame_count = 1;
  double r_cam = 75.0;
  double theta_cam = 1.5865;
  double phi_cam = 0.0;
  double observer_time = 0.0;
  double frame_dt = 8.0;
  double fov = 1.15;
  double lambda_max = 4000.0;
  double r_escape = 80.0;
  double exposure = 2.5;
  const char *output_path = "Output/accretion_disk.ppm";
  const char *hdr_output_path = "Output/accretion_disk_linear.pfm";
  const char *plot_output_path = "Output/accretion_disk_data.csv";
  const char *metadata_output_path = "Output/accretion_disk_metadata.json";
};

struct FishboneMoncriefTorusConfig {
  double r_in = 0.0;
  double r_center = 13.0;
  double gamma = 4.0 / 3.0;
  double rho_center = 1.0;
  double beta_inv = 0.03;
  double electron_density_scale = 4000.0;
  double thetae_scale = 120.0;
  double magnetic_scale = 18.0;
  double emissivity_scale = 2500.0;
  double opacity_scale = 40.0;
  double rho_cutoff = 1.0e-6;
};

struct FishboneMoncriefTorusModel {
  FishboneMoncriefTorusConfig config;
  double ell = 0.0;
  double Win = 0.0;
  double Kpoly = 0.0;
  double r_horizon = 0.0;
  double risco = 0.0;
};

struct ShadowFrameOutputs {
  std::string ppm;
  std::string pfm;
  std::string csv;
  std::string json;
};

struct SampleOffset {
  double x;
  double y;
};

constexpr SampleOffset kCenteredSample = {0.5, 0.5};
constexpr SampleOffset kSubpixelPattern[4] = {
    {0.25, 0.25},
    {0.75, 0.25},
    {0.25, 0.75},
    {0.75, 0.75},
};

static int read_env_int(const char *name, int fallback) {
  const char *value = std::getenv(name);
  if (!value || *value == '\0')
    return fallback;
  char *end = nullptr;
  long parsed = std::strtol(value, &end, 10);
  if (!end || *end != '\0')
    return fallback;
  if (parsed <= 0 || parsed > std::numeric_limits<int>::max())
    return fallback;
  return static_cast<int>(parsed);
}

static double read_env_double(const char *name, double fallback) {
  const char *value = std::getenv(name);
  if (!value || *value == '\0')
    return fallback;
  char *end = nullptr;
  double parsed = std::strtod(value, &end);
  if (!end || *end != '\0' || !std::isfinite(parsed))
    return fallback;
  return parsed;
}

static const char *read_env_path(const char *name, const char *fallback) {
  const char *value = std::getenv(name);
  if (!value || *value == '\0')
    return fallback;
  return value;
}

static ShadowRenderConfig load_shadow_render_config() {
  ShadowRenderConfig config;
  config.width = read_env_int("CURVATUREENGINE_SHADOW_WIDTH", config.width);
  config.height = read_env_int("CURVATUREENGINE_SHADOW_HEIGHT", config.height);
  config.samples_per_pixel =
      std::min(4, read_env_int("CURVATUREENGINE_SHADOW_SPP",
                               config.samples_per_pixel));
  config.frame_count =
      read_env_int("CURVATUREENGINE_SHADOW_FRAMES", config.frame_count);
  config.r_cam = read_env_double("CURVATUREENGINE_SHADOW_R_CAM", config.r_cam);
  config.theta_cam =
      read_env_double("CURVATUREENGINE_SHADOW_THETA_CAM", config.theta_cam);
  config.phi_cam =
      read_env_double("CURVATUREENGINE_SHADOW_PHI_CAM", config.phi_cam);
  config.observer_time =
      read_env_double("CURVATUREENGINE_SHADOW_TIME", config.observer_time);
  config.frame_dt =
      read_env_double("CURVATUREENGINE_SHADOW_FRAME_DT", config.frame_dt);
  config.fov = read_env_double("CURVATUREENGINE_SHADOW_FOV", config.fov);
  config.lambda_max =
      read_env_double("CURVATUREENGINE_SHADOW_LAMBDA_MAX", config.lambda_max);
  config.r_escape =
      read_env_double("CURVATUREENGINE_SHADOW_R_ESCAPE", config.r_escape);
  config.exposure =
      read_env_double("CURVATUREENGINE_SHADOW_EXPOSURE", config.exposure);
  config.output_path =
      read_env_path("CURVATUREENGINE_SHADOW_OUTPUT_PPM", config.output_path);
  config.hdr_output_path = read_env_path("CURVATUREENGINE_SHADOW_OUTPUT_PFM",
                                         config.hdr_output_path);
  config.plot_output_path = read_env_path("CURVATUREENGINE_SHADOW_OUTPUT_CSV",
                                          config.plot_output_path);
  config.metadata_output_path =
      read_env_path("CURVATUREENGINE_SHADOW_OUTPUT_JSON",
                    config.metadata_output_path);
  return config;
}

static FishboneMoncriefTorusConfig load_fm_torus_config() {
  FishboneMoncriefTorusConfig config;
  config.r_in = read_env_double("CURVATUREENGINE_TORUS_R_IN", config.r_in);
  config.r_center =
      read_env_double("CURVATUREENGINE_TORUS_R_CENTER", config.r_center);
  config.gamma = read_env_double("CURVATUREENGINE_TORUS_GAMMA", config.gamma);
  config.rho_center =
      read_env_double("CURVATUREENGINE_TORUS_RHO_CENTER", config.rho_center);
  config.beta_inv =
      read_env_double("CURVATUREENGINE_TORUS_BETA_INV", config.beta_inv);
  config.electron_density_scale = read_env_double(
      "CURVATUREENGINE_TORUS_NE_SCALE", config.electron_density_scale);
  config.thetae_scale = read_env_double("CURVATUREENGINE_TORUS_THETAE_SCALE",
                                        config.thetae_scale);
  config.magnetic_scale = read_env_double("CURVATUREENGINE_TORUS_B_SCALE",
                                          config.magnetic_scale);
  config.emissivity_scale = read_env_double(
      "CURVATUREENGINE_TORUS_J_SCALE", config.emissivity_scale);
  config.opacity_scale = read_env_double("CURVATUREENGINE_TORUS_A_SCALE",
                                         config.opacity_scale);
  config.rho_cutoff = read_env_double("CURVATUREENGINE_TORUS_RHO_CUTOFF",
                                      config.rho_cutoff);
  return config;
}

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

static inline double safe_pos(double x) { return (x > 1e-30) ? x : 1e-30; }

static inline double clamp01(double x) {
  if (x < 0.0)
    return 0.0;
  if (x > 1.0)
    return 1.0;
  return x;
}

static inline double wrap_angle(double x) {
  return std::atan2(std::sin(x), std::cos(x));
}

static inline void tone_map_rgb(double Ir, double Ig, double Ib, double &out_r,
                                double &out_g, double &out_b) {
  auto TM = [](double x) {
    x = std::max(0.0, x);
    double y = x / (1.0 + x);
    y = std::pow(y, 1.0 / 2.2);
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

static inline double compute_prograde_keplerian_omega(double r, double a_d) {
  return 1.0 / (std::pow(r, 1.5) + a_d);
}

static double compute_disk_dynamics_factor(
    double time_coord, double r, double theta, double phi,
    const FishboneMoncriefTorusModel &torus_model) {
  const double cylindrical_r = r * std::sin(theta);
  const double z = r * std::cos(theta);
  const double r_center = torus_model.config.r_center;
  const double hotspot_r = 0.82 * r_center;
  const double hotspot_sigma_r = 0.10 * r_center;
  const double hotspot_sigma_z = 0.020 * r_center;
  const double hotspot_phase_sigma = 0.35;
  const double hotspot_omega =
      compute_prograde_keplerian_omega(hotspot_r, static_cast<double>(a));
  const double hotspot_phase = wrap_angle(phi - hotspot_omega * time_coord);
  const double hotspot =
      std::exp(-0.5 * std::pow((cylindrical_r - hotspot_r) /
                                   std::max(hotspot_sigma_r, 1.0e-6),
                               2.0) -
               0.5 * std::pow(z / std::max(hotspot_sigma_z, 1.0e-6), 2.0) -
               0.5 * std::pow(hotspot_phase / hotspot_phase_sigma, 2.0));

  const double pattern_omega = 0.55 * hotspot_omega;
  const double spiral_phase =
      2.0 * (phi - pattern_omega * time_coord) -
      5.0 * std::log(std::max(cylindrical_r, 1.0e-6) /
                     std::max(r_center, 1.0e-6));
  const double spiral =
      0.5 * (1.0 + std::cos(spiral_phase)) *
      std::exp(-0.5 * std::pow((cylindrical_r - r_center) /
                                   std::max(0.45 * r_center, 1.0e-6),
                               2.0));

  const double local_omega = compute_prograde_keplerian_omega(
      std::max(cylindrical_r, torus_model.r_horizon + 0.5),
      static_cast<double>(a));
  const double turbulence =
      0.5 * (std::sin(3.0 * phi - 1.15 * local_omega * time_coord +
                      0.7 * std::log(std::max(cylindrical_r, 1.0e-6))) +
             std::cos(5.0 * phi - 0.75 * hotspot_omega * time_coord));

  const double factor =
      1.0 + 1.5 * hotspot + 0.45 * spiral + 0.12 * turbulence;
  return std::clamp(factor, 0.25, 3.25);
}

static inline bool compute_fluid_state_from_omega(double g_tt, double g_tph,
                                                  double g_phph, double omega,
                                                  double &ut, double &uph,
                                                  double &u_cov_t,
                                                  double &u_cov_ph) {
  double norm_u = -(g_tt + 2.0 * g_tph * omega + g_phph * omega * omega);
  if (!(norm_u > 0.0))
    return false;
  ut = 1.0 / std::sqrt(norm_u);
  uph = omega * ut;
  u_cov_t = g_tt * ut + g_tph * uph;
  u_cov_ph = g_tph * ut + g_phph * uph;
  return true;
}

static bool sample_metric_components(double r, double theta, double &g_tt,
                                     double &g_tph, double &g_phph) {
  alignas(32) double Xt[4] = {0.0, 0.0, 0.0, 0.0};
  alignas(32) double Xr[4] = {r, r, r, r};
  alignas(32) double Xth[4] = {theta, theta, theta, theta};
  alignas(32) double Xph[4] = {0.0, 0.0, 0.0, 0.0};

  VEC_TYPE x_mid[NDIM] = {load(Xt), load(Xr), load(Xth), load(Xph)};
  VEC_TYPE g[NDIM][NDIM], g_inv[NDIM][NDIM];
  MetricSIMD::calculate_metric_avx(x_mid, g, g_inv);

  alignas(32) double lane[4];
  store(lane, g[0][0]);
  g_tt = lane[0];
  store(lane, g[0][3]);
  g_tph = lane[0];
  store(lane, g[3][3]);
  g_phph = lane[0];
  return std::isfinite(g_tt) && std::isfinite(g_tph) && std::isfinite(g_phph);
}

static inline bool compute_constant_ell_fluid_from_metric(
    double g_tt, double g_tph, double g_phph, double ell, double &W,
    double &ut, double &uph, double &u_cov_t, double &u_cov_ph);

static bool sample_constant_ell_potential_equatorial(double r, double ell,
                                                     double &W) {
  double g_tt = 0.0;
  double g_tph = 0.0;
  double g_phph = 0.0;
  double ut = 0.0;
  double uph = 0.0;
  double u_cov_t = 0.0;
  double u_cov_ph = 0.0;
  if (!sample_metric_components(r, kPi / 2.0, g_tt, g_tph, g_phph)) {
    return false;
  }
  return compute_constant_ell_fluid_from_metric(g_tt, g_tph, g_phph, ell, W,
                                                ut, uph, u_cov_t, u_cov_ph);
}

static bool find_bound_surface_radius(double ell, double r_min, double r_max,
                                      double &surface_radius) {
  constexpr int kScanSteps = 512;
  double previous_r = 0.0;
  double previous_W = 0.0;
  bool have_previous = false;

  for (int step = 0; step <= kScanSteps; ++step) {
    const double t = static_cast<double>(step) / static_cast<double>(kScanSteps);
    const double r = r_min + (r_max - r_min) * t;
    double W = 0.0;
    if (!sample_constant_ell_potential_equatorial(r, ell, W))
      continue;

    if (have_previous && previous_W > 0.0 && W <= 0.0) {
      double lo = previous_r;
      double hi = r;
      for (int iter = 0; iter < 48; ++iter) {
        const double mid = 0.5 * (lo + hi);
        double W_mid = 0.0;
        if (!sample_constant_ell_potential_equatorial(mid, ell, W_mid))
          return false;
        if (W_mid > 0.0) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      surface_radius = 0.5 * (lo + hi);
      return true;
    }

    previous_r = r;
    previous_W = W;
    have_previous = true;
  }

  return false;
}

static inline bool compute_constant_ell_fluid_from_metric(
    double g_tt, double g_tph, double g_phph, double ell, double &W,
    double &ut, double &uph, double &u_cov_t, double &u_cov_ph) {
  const double denom = g_phph + ell * g_tph;
  if (std::abs(denom) < 1e-30)
    return false;

  const double omega = (-g_tph - ell * g_tt) / denom;
  if (!compute_fluid_state_from_omega(g_tt, g_tph, g_phph, omega, ut, uph,
                                      u_cov_t, u_cov_ph)) {
    return false;
  }

  const double minus_ut = -u_cov_t;
  if (!(minus_ut > 0.0))
    return false;
  W = std::log(minus_ut);
  return std::isfinite(W);
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

static bool build_fm_torus_model(double spin,
                                 const FishboneMoncriefTorusConfig &input,
                                 FishboneMoncriefTorusModel &model) {
  model = {};
  model.config = input;
  model.r_horizon = 1.0 + std::sqrt(std::max(0.0, 1.0 - spin * spin));
  model.risco = compute_risco(spin);
  const bool user_defined_r_in = input.r_in > 0.0;

  if (!(model.config.gamma > 1.0) || !(model.config.rho_center > 0.0) ||
      !(model.config.beta_inv >= 0.0) ||
      !(model.config.electron_density_scale >= 0.0) ||
      !(model.config.thetae_scale >= 0.0) ||
      !(model.config.magnetic_scale >= 0.0) ||
      !(model.config.emissivity_scale >= 0.0) ||
      !(model.config.opacity_scale >= 0.0) ||
      !(model.config.rho_cutoff >= 0.0)) {
    return false;
  }

  if (!(model.config.r_in > 0.0)) {
    model.config.r_in = std::max(model.r_horizon + 0.4, model.risco + 0.6);
  } else {
    model.config.r_in = std::max(model.config.r_in, model.r_horizon + 1e-3);
  }

  if (!(model.config.r_center > 0.0)) {
    model.config.r_center =
        std::max(model.config.r_in + 8.0, model.risco + 10.0);
  } else if (model.config.r_center <= model.config.r_in + 0.5) {
    model.config.r_center = model.config.r_in + 2.0;
  }

  double g_tt_center, g_tph_center, g_phph_center;
  if (!sample_metric_components(model.config.r_center, kPi / 2.0, g_tt_center,
                                g_tph_center, g_phph_center)) {
    return false;
  }

  const double omega_center =
      compute_prograde_keplerian_omega(model.config.r_center, spin);
  double ut_center_orbit, uph_center_orbit, u_cov_t_center_orbit,
      u_cov_ph_center_orbit;
  if (!compute_fluid_state_from_omega(g_tt_center, g_tph_center,
                                      g_phph_center, omega_center,
                                      ut_center_orbit, uph_center_orbit,
                                      u_cov_t_center_orbit,
                                      u_cov_ph_center_orbit)) {
    return false;
  }
  (void)ut_center_orbit;
  (void)uph_center_orbit;
  model.ell = -u_cov_ph_center_orbit / u_cov_t_center_orbit;

  double W_center_equatorial = 0.0;
  if (!sample_constant_ell_potential_equatorial(model.config.r_center,
                                                model.ell,
                                                W_center_equatorial)) {
    return false;
  }
  if (!(W_center_equatorial < 0.0) || !std::isfinite(W_center_equatorial))
    return false;

  double r_bound_surface = 0.0;
  if (!find_bound_surface_radius(model.ell, model.r_horizon + 0.02,
                                 model.config.r_center - 1.0e-3,
                                 r_bound_surface)) {
    return false;
  }

  // Shift the default torus farther out to widen the disk radially, while
  // keeping a relatively compact cross-section so it does not puff up.
  const double compact_offset =
      0.25 * (model.config.r_center - r_bound_surface);
  if (!user_defined_r_in) {
    model.config.r_in = r_bound_surface + compact_offset;
  } else if (model.config.r_in <= r_bound_surface) {
    model.config.r_in = r_bound_surface + 0.05 *
                                              (model.config.r_center -
                                               r_bound_surface);
  }

  double g_tt_inner, g_tph_inner, g_phph_inner;
  if (!sample_metric_components(model.config.r_in, kPi / 2.0, g_tt_inner,
                                g_tph_inner, g_phph_inner)) {
    return false;
  }

  double W_inner, W_center, ut_dummy, uph_dummy, u_cov_t_dummy,
      u_cov_ph_dummy;
  if (!compute_constant_ell_fluid_from_metric(
          g_tt_inner, g_tph_inner, g_phph_inner, model.ell, W_inner, ut_dummy,
          uph_dummy, u_cov_t_dummy, u_cov_ph_dummy)) {
    return false;
  }
  if (!compute_constant_ell_fluid_from_metric(g_tt_center, g_tph_center,
                                              g_phph_center, model.ell,
                                              W_center, ut_dummy, uph_dummy,
                                              u_cov_t_dummy, u_cov_ph_dummy)) {
    return false;
  }
  (void)ut_dummy;
  (void)uph_dummy;

  const double delta_W = W_inner - W_center;
  if (!(W_inner < 0.0) || !(W_inner > W_center_equatorial) ||
      !(delta_W > 0.0) || !std::isfinite(delta_W))
    return false;

  const double h_center = std::exp(std::min(delta_W, 60.0));
  if (!(h_center > 1.0))
    return false;

  const double rho_center_pow =
      std::pow(model.config.rho_center, model.config.gamma - 1.0);
  const double denom =
      model.config.gamma * std::max(rho_center_pow, 1.0e-30);
  model.Kpoly =
      (h_center - 1.0) * (model.config.gamma - 1.0) / std::max(denom, 1.0e-30);
  model.Win = W_inner;

  return std::isfinite(model.ell) && std::isfinite(model.Win) &&
         std::isfinite(model.Kpoly) && model.Kpoly > 0.0;
}

static inline bool fm_torus_state_from_metric(
    double g_tt, double g_tph, double g_phph, double ell, double Win,
    double Gamma, double Kpoly, double &rho, double &p, double &ut,
    double &uph, double &u_cov_t, double &u_cov_ph) {
  double W = 0.0;
  if (!compute_constant_ell_fluid_from_metric(g_tt, g_tph, g_phph, ell, W, ut,
                                              uph, u_cov_t, u_cov_ph)) {
    return false;
  }

  const double delta_W = Win - W;
  if (!(delta_W > 0.0))
    return false;

  const double h = std::exp(std::min(delta_W, 60.0));
  if (!(h > 1.0)) {
    rho = 0.0;
    p = 0.0;
    return false;
  }

  const double x = (h - 1.0) * (Gamma - 1.0) / safe_pos(Gamma * Kpoly);
  if (!(x > 0.0))
    return false;

  rho = std::pow(x, 1.0 / (Gamma - 1.0));
  p = Kpoly * std::pow(rho, Gamma);
  return std::isfinite(rho) && std::isfinite(p) && rho > 0.0 && p > 0.0;
}

static inline void accumulate_grrt_step(
    const VEC_TYPE x_prev[NDIM], const VEC_TYPE x_curr[NDIM],
    const VEC_TYPE k_curr[NDIM], const double nu_obs[4], double dlambda,
    const FishboneMoncriefTorusModel &torus_model, double Ibar[4],
    double tau[4], double WTe[4], double Wsum[4]) {
  constexpr double kThetaeToKelvin = 5.930000000e9;
  constexpr double kEmissionHeightFactor = 0.035;
  constexpr double kRadialSigmaFactor = 0.32;

  alignas(32) double tP[4], tC[4], rP[4], rC[4], thP[4], thC[4], phP[4], phC[4];
  store(tP, x_prev[0]);
  store(tC, x_curr[0]);
  store(rP, x_prev[1]);
  store(rC, x_curr[1]);
  store(thP, x_prev[2]);
  store(thC, x_curr[2]);
  store(phP, x_prev[3]);
  store(phC, x_curr[3]);

  alignas(32) double kt[4], kr[4], kth[4], kph[4];
  store(kt, k_curr[0]);
  store(kr, k_curr[1]);
  store(kth, k_curr[2]);
  store(kph, k_curr[3]);

  alignas(32) double tM[4], rM[4], thM[4], phM[4];
  for (int i = 0; i < 4; ++i) {
    tM[i] = 0.5 * (tP[i] + tC[i]);
    rM[i] = 0.5 * (rP[i] + rC[i]);
    thM[i] = 0.5 * (thP[i] + thC[i]);
    phM[i] = 0.5 * (phP[i] + phC[i]);
  }

  alignas(32) double Xt[4] = {0}, Xph[4] = {0};
  alignas(32) double Xr[4] = {rM[0], rM[1], rM[2], rM[3]};
  alignas(32) double Xth[4] = {thM[0], thM[1], thM[2], thM[3]};

  VEC_TYPE X_mid[NDIM] = {load(Xt), load(Xr), load(Xth), load(Xph)};
  VEC_TYPE gS[NDIM][NDIM], ginvS[NDIM][NDIM];
  MetricSIMD::calculate_metric_avx(X_mid, gS, ginvS);

  alignas(32) double gtt[4], grr[4], gthth[4], gtph[4], gphph[4];
  store(gtt, gS[0][0]);
  store(grr, gS[1][1]);
  store(gthth, gS[2][2]);
  store(gtph, gS[0][3]);
  store(gphph, gS[3][3]);

  for (int i = 0; i < 4; ++i) {
    const double r = rM[i];
    const double theta = thM[i];
    const double phi = phM[i];
    const double time_coord = tM[i];
    if (r <= torus_model.r_horizon + 1e-3)
      continue;

    double rho = 0.0;
    double p = 0.0;
    double ut = 0.0;
    double uph = 0.0;
    double u_cov_t = 0.0;
    double u_cov_ph = 0.0;
    if (!fm_torus_state_from_metric(gtt[i], gtph[i], gphph[i], torus_model.ell,
                                    torus_model.Win, torus_model.config.gamma,
                                    torus_model.Kpoly, rho, p, ut, uph,
                                    u_cov_t, u_cov_ph)) {
      continue;
    }
    (void)ut;
    (void)uph;

    if (rho < torus_model.config.rho_cutoff || p <= 0.0)
      continue;

    const double nu_em =
        -(u_cov_t * kt[i] + u_cov_ph * kph[i]);
    if (!(nu_em > 1e-12))
      continue;

    const double redshift = nu_obs[i] / nu_em;
    const double cylindrical_r = r * std::sin(theta);
    const double z = r * std::cos(theta);
    const double emission_height =
        kEmissionHeightFactor * torus_model.config.r_center;
    const double vertical_taper =
        std::exp(-0.5 * (z / std::max(emission_height, 1.0e-6)) *
                 (z / std::max(emission_height, 1.0e-6)));
    if (vertical_taper < 1.0e-5)
      continue;

    const double radial_sigma =
        kRadialSigmaFactor * torus_model.config.r_center;
    const double radial_weight =
        0.35 +
        0.65 * std::exp(-0.5 * std::pow(
                                   (cylindrical_r - torus_model.config.r_center) /
                                       std::max(radial_sigma, 1.0e-6),
                                   2.0));
    const double dynamics_factor =
        compute_disk_dynamics_factor(time_coord, r, theta, phi, torus_model);

    double g_cov[4][4] = {
        {gtt[i], 0.0, 0.0, gtph[i]},
        {0.0, grr[i], 0.0, 0.0},
        {0.0, 0.0, gthth[i], 0.0},
        {gtph[i], 0.0, 0.0, gphph[i]},
    };
    double u_vec[4] = {ut, 0.0, 0.0, uph};
    double photon_vec[4] = {kt[i], kr[i], kth[i], kph[i]};
    double toroidal_b[4] = {0.0, 0.0, 0.0, 1.0};
    project_orthogonal_to_u(g_cov, u_vec, toroidal_b);
    normalize_spacelike(g_cov, toroidal_b);

    double photon_dir[4];
    for (int mu = 0; mu < 4; ++mu)
      photon_dir[mu] = photon_vec[mu] / nu_em - u_vec[mu];
    project_orthogonal_to_u(g_cov, u_vec, photon_dir);
    normalize_spacelike(g_cov, photon_dir);

    const double cos_pitch = std::clamp(dot_g(g_cov, photon_dir, toroidal_b),
                                        -1.0, 1.0);
    const double sin_pitch =
        std::sqrt(std::max(0.0, 1.0 - cos_pitch * cos_pitch));
    const double synchrotron_anisotropy = 0.2 + 0.8 * sin_pitch;
    const double doppler_boost =
        std::clamp(std::pow(redshift, 0.85), 0.35, 3.0);

    const double ne =
        torus_model.config.electron_density_scale * rho * dynamics_factor;
    const double thetae = std::max(
        1.0e-4,
        torus_model.config.thetae_scale * p / safe_pos(rho) *
            (0.9 + 0.1 * dynamics_factor));
    const double B =
        torus_model.config.magnetic_scale *
        std::sqrt(std::max(0.0, torus_model.config.beta_inv * p)) *
        std::sqrt(dynamics_factor);
    if (!(ne > 0.0) || !(B > 0.0))
      continue;

    double jnu = 0.0;
    double anu = 0.0;
    emiss_abs_synch_approx(ne, thetae, B, nu_em, jnu, anu);
    const double emissivity_modifier =
        vertical_taper * radial_weight * synchrotron_anisotropy *
        doppler_boost;
    const double opacity_modifier =
        vertical_taper * radial_weight * std::sqrt(synchrotron_anisotropy);
    jnu *= torus_model.config.emissivity_scale * emissivity_modifier;
    anu *= torus_model.config.opacity_scale * opacity_modifier;

    if (!(jnu > 0.0) && !(anu > 0.0))
      continue;

    const double Jinv = jnu / (nu_em * nu_em);
    const double Ainv = anu * nu_em;

    const double I0 = Ibar[i];
    const double tau0 = tau[i];

    rt_step_exact(dlambda, Jinv, Ainv, Ibar[i], tau[i]);

    const double dI = Ibar[i] - I0;
    if (dI > 0.0) {
      const double w = dI * std::exp(-tau0);
      const double Te_obs = thetae * kThetaeToKelvin * redshift;
      WTe[i] += Te_obs * w;
      Wsum[i] += w;
    }
  }
}

static PixelResult
geodesic_raytrace_physical(VEC_TYPE x[NDIM], VEC_TYPE v[NDIM],
                           const double uk_obs[4], float lambda_max,
                           ChristoffelEvalFnVec eval_fn, void *ctx,
                           double r_horizon, double r_escape,
                           const FishboneMoncriefTorusModel &torus_model,
                           double exposure) {
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

    accumulate_grrt_step(prev_x, x, v, nu_obs, step_val, torus_model, Ibar,
                         tau, WTe, Wsum);
  }

  PixelResult res;
  for (int i = 0; i < 4; ++i) {
    const bool has_emission = Wsum[i] > 0.0;
    double Te_eff = has_emission ? (WTe[i] / Wsum[i]) : 0.0;
    double Iobs = Ibar[i] * nu_obs[i] * nu_obs[i] * nu_obs[i];
    const double color_temperature =
        has_emission
            ? std::clamp(3500.0 + 4500.0 * std::log10(1.0 + Te_eff / 1.0e6),
                         1000.0, 40000.0)
            : 9000.0;
    double cr, cg, cb;
    blackbody_rgb(color_temperature, cr, cg, cb);
    double Ir = Iobs * cr * exposure;
    double Ig = Iobs * cg * exposure;
    double Ib = Iobs * cb * exposure;

    res.linear_r[i] = Ir;
    res.linear_g[i] = Ig;
    res.linear_b[i] = Ib;
    res.intensity[i] = Iobs;
    res.optical_depth[i] = tau[i];
    res.effective_temperature[i] = Te_eff;
  }
  return res;
}

bool ensure_parent_directory(const char *filename) {
  if (!filename || *filename == '\0')
    return false;
  const std::filesystem::path path(filename);
  if (!path.has_parent_path())
    return true;
  std::error_code ec;
  std::filesystem::create_directories(path.parent_path(), ec);
  return !ec;
}

std::string format_frame_output_path(const char *filename, int frame_index,
                                     int frame_count) {
  if (!filename)
    return "";
  if (frame_count <= 1)
    return filename;

  const std::filesystem::path path(filename);
  const std::string stem = path.stem().string();
  const std::string extension = path.extension().string();

  char suffix[16];
  std::snprintf(suffix, sizeof(suffix), "_%04d", frame_index);

  const std::filesystem::path numbered =
      path.parent_path() / (stem + suffix + extension);
  return numbered.string();
}

ShadowFrameOutputs build_frame_outputs(const ShadowRenderConfig &config,
                                       int frame_index) {
  ShadowFrameOutputs outputs;
  outputs.ppm =
      format_frame_output_path(config.output_path, frame_index, config.frame_count);
  outputs.pfm = format_frame_output_path(config.hdr_output_path, frame_index,
                                         config.frame_count);
  outputs.csv = format_frame_output_path(config.plot_output_path, frame_index,
                                         config.frame_count);
  outputs.json = format_frame_output_path(config.metadata_output_path,
                                          frame_index, config.frame_count);
  return outputs;
}

std::string build_sequence_pattern(const char *filename) {
  if (!filename)
    return "";

  const std::filesystem::path path(filename);
  const std::string stem = path.stem().string();
  const std::string extension = path.extension().string();
  return (path.parent_path() / (stem + "_%04d" + extension)).string();
}

std::string build_sequence_video_output(const char *filename) {
  if (!filename)
    return "Output/accretion_disk.mp4";

  const std::filesystem::path path(filename);
  return (path.parent_path() / (path.stem().string() + ".mp4")).string();
}

bool write_ppm_image(const char *filename, int width, int height,
                     const std::vector<uint8_t> &rgb_buffer) {
  if (!ensure_parent_directory(filename))
    return false;
  FILE *fp = fopen(filename, "wb");
  if (!fp)
    return false;
  fprintf(fp, "P6\n%d %d\n255\n", width, height);
  fwrite(rgb_buffer.data(), 1, rgb_buffer.size(), fp);
  fclose(fp);
  return true;
}

bool write_pfm_image(const char *filename, int width, int height,
                     const std::vector<double> &rgb_buffer) {
  if (!ensure_parent_directory(filename))
    return false;
  FILE *fp = fopen(filename, "wb");
  if (!fp)
    return false;

  fprintf(fp, "PF\n%d %d\n-1.0\n", width, height);
  std::vector<float> row(static_cast<std::size_t>(width) * 3U);
  for (int y = height - 1; y >= 0; --y) {
    for (int x = 0; x < width; ++x) {
      const std::size_t src = (static_cast<std::size_t>(y) * width + x) * 3U;
      const std::size_t dst = static_cast<std::size_t>(x) * 3U;
      row[dst + 0] = static_cast<float>(rgb_buffer[src + 0]);
      row[dst + 1] = static_cast<float>(rgb_buffer[src + 1]);
      row[dst + 2] = static_cast<float>(rgb_buffer[src + 2]);
    }
    fwrite(row.data(), sizeof(float), row.size(), fp);
  }
  fclose(fp);
  return true;
}

bool write_shadow_plot_data_csv(const char *filename,
                                const ShadowRenderConfig &config,
                                const std::vector<double> &rgb_buffer,
                                const std::vector<double> &intensity_buffer,
                                const std::vector<double> &tau_buffer,
                                const std::vector<double> &temperature_buffer) {
  if (!ensure_parent_directory(filename))
    return false;

  std::ofstream file(filename);
  if (!file.is_open())
    return false;

  const double aspect =
      static_cast<double>(config.width) / static_cast<double>(config.height);
  file << std::scientific << std::setprecision(12);
  file << "pixel_x,pixel_y,screen_x,screen_y,linear_r,linear_g,linear_b,"
          "intensity,optical_depth,effective_temperature\n";

  for (int y = 0; y < config.height; ++y) {
    const double screen_y =
        -(2.0 * ((y + 0.5) / static_cast<double>(config.height)) - 1.0);
    for (int x = 0; x < config.width; ++x) {
      const double screen_x =
          (2.0 * ((x + 0.5) / static_cast<double>(config.width)) - 1.0) *
          aspect;
      const std::size_t pixel_idx =
          static_cast<std::size_t>(y) * config.width + x;
      const std::size_t rgb_idx = pixel_idx * 3U;
      file << x << ',' << y << ',' << screen_x << ',' << screen_y << ','
           << rgb_buffer[rgb_idx + 0] << ',' << rgb_buffer[rgb_idx + 1] << ','
           << rgb_buffer[rgb_idx + 2] << ',' << intensity_buffer[pixel_idx]
           << ',' << tau_buffer[pixel_idx] << ','
           << temperature_buffer[pixel_idx] << '\n';
    }
  }

  return true;
}

bool write_shadow_metadata_json(const char *filename,
                                const ShadowRenderConfig &config,
                                const FishboneMoncriefTorusModel &torus_model,
                                int frame_index) {
  if (!ensure_parent_directory(filename))
    return false;

  std::ofstream file(filename);
  if (!file.is_open())
    return false;

  file << std::fixed << std::setprecision(6);
  file << "{\n";
  file << "  \"model\": \"fishbone_moncrief\",\n";
  file << "  \"width\": " << config.width << ",\n";
  file << "  \"height\": " << config.height << ",\n";
  file << "  \"samples_per_pixel\": " << config.samples_per_pixel << ",\n";
  file << "  \"frame_index\": " << frame_index << ",\n";
  file << "  \"frame_count\": " << config.frame_count << ",\n";
  file << "  \"observer_time\": " << config.observer_time << ",\n";
  file << "  \"frame_dt\": " << config.frame_dt << ",\n";
  file << "  \"spin\": " << static_cast<double>(a) << ",\n";
  file << "  \"r_cam\": " << config.r_cam << ",\n";
  file << "  \"theta_cam\": " << config.theta_cam << ",\n";
  file << "  \"phi_cam\": " << config.phi_cam << ",\n";
  file << "  \"fov\": " << config.fov << ",\n";
  file << "  \"lambda_max\": " << config.lambda_max << ",\n";
  file << "  \"r_escape\": " << config.r_escape << ",\n";
  file << "  \"exposure\": " << config.exposure << ",\n";
  file << "  \"ppm_output\": \"" << config.output_path << "\",\n";
  file << "  \"pfm_output\": \"" << config.hdr_output_path << "\",\n";
  file << "  \"csv_output\": \"" << config.plot_output_path << "\",\n";
  file << "  \"torus\": {\n";
  file << "    \"r_horizon\": " << torus_model.r_horizon << ",\n";
  file << "    \"risco\": " << torus_model.risco << ",\n";
  file << "    \"r_in\": " << torus_model.config.r_in << ",\n";
  file << "    \"r_center\": " << torus_model.config.r_center << ",\n";
  file << "    \"gamma\": " << torus_model.config.gamma << ",\n";
  file << "    \"rho_center\": " << torus_model.config.rho_center << ",\n";
  file << "    \"beta_inv\": " << torus_model.config.beta_inv << ",\n";
  file << "    \"electron_density_scale\": "
       << torus_model.config.electron_density_scale << ",\n";
  file << "    \"thetae_scale\": " << torus_model.config.thetae_scale
       << ",\n";
  file << "    \"magnetic_scale\": " << torus_model.config.magnetic_scale
       << ",\n";
  file << "    \"emissivity_scale\": " << torus_model.config.emissivity_scale
       << ",\n";
  file << "    \"opacity_scale\": " << torus_model.config.opacity_scale
       << ",\n";
  file << "    \"rho_cutoff\": " << torus_model.config.rho_cutoff << ",\n";
  file << "    \"ell\": " << torus_model.ell << ",\n";
  file << "    \"Win\": " << torus_model.Win << ",\n";
  file << "    \"Kpoly\": " << torus_model.Kpoly << "\n";
  file << "  }\n";
  file << "}\n";
  return true;
}

bool is_valid_render_config(const ShadowRenderConfig &config) {
  return config.width > 0 && config.height > 0 && config.r_cam > 0.0 &&
         config.lambda_max > 0.0 && config.r_escape > 0.0 &&
         config.samples_per_pixel > 0 && config.frame_count > 0 &&
         config.frame_dt >= 0.0 && std::isfinite(config.observer_time) &&
         std::isfinite(config.frame_dt) && config.output_path != nullptr &&
         config.hdr_output_path != nullptr &&
         config.plot_output_path != nullptr &&
         config.metadata_output_path != nullptr;
}

bool is_valid_fm_torus_config(const FishboneMoncriefTorusConfig &config) {
  return config.gamma > 1.0 && config.rho_center > 0.0 &&
         config.beta_inv >= 0.0 && config.electron_density_scale >= 0.0 &&
         config.thetae_scale >= 0.0 && config.magnetic_scale >= 0.0 &&
         config.emissivity_scale >= 0.0 && config.opacity_scale >= 0.0 &&
         config.rho_cutoff >= 0.0;
}

void compute_observer_frequency(const Tetrad &tetrad, const double g_cov[4][4],
                                const double kt[4], const double kr[4],
                                const double kth[4], const double kph[4],
                                double uk_obs[4]) {
  double u_obs_cov[4] = {0.0, 0.0, 0.0, 0.0};
  for (int mu = 0; mu < 4; ++mu) {
    double sum = 0.0;
    for (int nu = 0; nu < 4; ++nu) {
      sum += g_cov[mu][nu] * tetrad.e[0][nu];
    }
    u_obs_cov[mu] = sum;
  }

  for (int i = 0; i < 4; ++i) {
    uk_obs[i] = u_obs_cov[0] * kt[i] + u_obs_cov[1] * kr[i] +
                u_obs_cov[2] * kth[i] + u_obs_cov[3] * kph[i];
  }
}

PixelResult render_packet(int y, int packet_index,
                          const ShadowRenderConfig &config, double aspect,
                          double sample_x, double sample_y,
                          const FishboneMoncriefTorusModel &torus_model) {
  const int x_base = packet_index * 4;
  const double sy =
      -(2.0 * ((y + sample_y) / static_cast<double>(config.height)) - 1.0);

  alignas(32) double sx4[4];
  for (int lane = 0; lane < 4; ++lane) {
    const int x = std::min(x_base + lane, config.width - 1);
    const double sx =
        2.0 * ((x + sample_x) / static_cast<double>(config.width)) - 1.0;
    sx4[lane] = sx * aspect;
  }

  VEC_TYPE x_cam[NDIM] = {broadcast(config.observer_time),
                          broadcast(config.r_cam), broadcast(config.theta_cam),
                          broadcast(config.phi_cam)};
  VEC_TYPE g_simd[NDIM][NDIM], ginv_simd[NDIM][NDIM];
  MetricSIMD::calculate_metric_avx(x_cam, g_simd, ginv_simd);

  alignas(32) double cov_lane[4], inv_lane[4];
  double g_cov[4][4], g_inv[4][4];
  for (int mu = 0; mu < 4; ++mu) {
    for (int nu = 0; nu < 4; ++nu) {
      store(cov_lane, g_simd[mu][nu]);
      store(inv_lane, ginv_simd[mu][nu]);
      g_cov[mu][nu] = cov_lane[0];
      g_inv[mu][nu] = inv_lane[0];
    }
  }

  Tetrad tetrad;
  build_camera_tetrad_zamo(g_cov, g_inv, tetrad);

  const double sy4[4] = {sy, sy, sy, sy};
  double k_out[4][4];
  make_camera_rays_4(tetrad, sx4, sy4, config.fov, k_out);

  alignas(32) double kt[4], kr[4], kth[4], kph[4];
  for (int i = 0; i < 4; ++i) {
    kt[i] = k_out[i][0];
    kr[i] = k_out[i][1];
    kth[i] = k_out[i][2];
    kph[i] = k_out[i][3];
  }

  double uk_obs[4];
  compute_observer_frequency(tetrad, g_cov, kt, kr, kth, kph, uk_obs);

  const double r_stop = torus_model.r_horizon + 0.05;
  VEC_TYPE v_avx[NDIM] = {load(kt), load(kr), load(kth), load(kph)};
  return geodesic_raytrace_physical(x_cam, v_avx, uk_obs,
                                    static_cast<float>(config.lambda_max),
                                    evaluate_christoffel_native_avx, nullptr,
                                    r_stop, config.r_escape, torus_model,
                                    config.exposure);
}

int render_shadow_frame(const ShadowRenderConfig &base_config,
                        const FishboneMoncriefTorusModel &torus_model,
                        int frame_index) {
  ShadowRenderConfig config = base_config;
  ShadowFrameOutputs outputs = build_frame_outputs(base_config, frame_index);
  config.observer_time =
      base_config.observer_time +
      static_cast<double>(frame_index) * base_config.frame_dt;
  config.output_path = outputs.ppm.c_str();
  config.hdr_output_path = outputs.pfm.c_str();
  config.plot_output_path = outputs.csv.c_str();
  config.metadata_output_path = outputs.json.c_str();

  std::vector<uint8_t> image_buffer(config.width * config.height * 3);
  std::vector<double> hdr_buffer(static_cast<std::size_t>(config.width) *
                                 config.height * 3U, 0.0);
  std::vector<double> intensity_buffer(static_cast<std::size_t>(config.width) *
                                       config.height, 0.0);
  std::vector<double> tau_buffer(static_cast<std::size_t>(config.width) *
                                 config.height, 0.0);
  std::vector<double> temperature_buffer(
      static_cast<std::size_t>(config.width) * config.height, 0.0);
  const double aspect =
      static_cast<double>(config.width) / static_cast<double>(config.height);
  const int packets_x = (config.width + 3) / 4;
  const int sample_count = std::max(1, std::min(config.samples_per_pixel, 4));

#pragma omp parallel for schedule(dynamic)
  for (int y = 0; y < config.height; ++y) {
    for (int p = 0; p < packets_x; ++p) {
      const int x_base = p * 4;
      PixelResult accum{};

      for (int sample = 0; sample < sample_count; ++sample) {
        const SampleOffset offset =
            (sample_count == 1) ? kCenteredSample : kSubpixelPattern[sample];
        const PixelResult res = render_packet(y, p, config, aspect, offset.x,
                                              offset.y, torus_model);
        for (int l = 0; l < 4; ++l) {
          accum.linear_r[l] += res.linear_r[l];
          accum.linear_g[l] += res.linear_g[l];
          accum.linear_b[l] += res.linear_b[l];
          accum.intensity[l] += res.intensity[l];
          accum.optical_depth[l] += res.optical_depth[l];
          accum.effective_temperature[l] += res.effective_temperature[l];
        }
      }

      const double inv_samples = 1.0 / static_cast<double>(sample_count);

      for (int l = 0; l < 4; ++l) {
        const int x = x_base + l;
        if (x >= config.width) {
          continue;
        }
        const std::size_t pixel_idx =
            static_cast<std::size_t>(y) * config.width + x;
        const int px_idx = (y * config.width + x) * 3;
        const std::size_t hdr_idx = pixel_idx * 3U;
        const double linear_r = accum.linear_r[l] * inv_samples;
        const double linear_g = accum.linear_g[l] * inv_samples;
        const double linear_b = accum.linear_b[l] * inv_samples;
        const double intensity = accum.intensity[l] * inv_samples;
        const double optical_depth = accum.optical_depth[l] * inv_samples;
        const double effective_temperature =
            accum.effective_temperature[l] * inv_samples;

        hdr_buffer[hdr_idx + 0] = linear_r;
        hdr_buffer[hdr_idx + 1] = linear_g;
        hdr_buffer[hdr_idx + 2] = linear_b;
        intensity_buffer[pixel_idx] = intensity;
        tau_buffer[pixel_idx] = optical_depth;
        temperature_buffer[pixel_idx] = effective_temperature;

        double display_r, display_g, display_b;
        tone_map_rgb(linear_r, linear_g, linear_b, display_r, display_g,
                     display_b);
        image_buffer[px_idx + 0] = static_cast<uint8_t>(display_r);
        image_buffer[px_idx + 1] = static_cast<uint8_t>(display_g);
        image_buffer[px_idx + 2] = static_cast<uint8_t>(display_b);
      }
    }
  }

  const bool ppm_ok =
      write_ppm_image(config.output_path, config.width, config.height,
                      image_buffer);
  const bool pfm_ok =
      write_pfm_image(config.hdr_output_path, config.width, config.height,
                      hdr_buffer);
  const bool csv_ok = write_shadow_plot_data_csv(
      config.plot_output_path, config, hdr_buffer, intensity_buffer, tau_buffer,
      temperature_buffer);
  const bool json_ok = write_shadow_metadata_json(
      config.metadata_output_path, config, torus_model, frame_index);

  if (!ppm_ok || !pfm_ok || !csv_ok || !json_ok) {
    std::fprintf(stderr,
                 "Failed to export shadow outputs (ppm=%d pfm=%d csv=%d json=%d).\n",
                 ppm_ok ? 1 : 0, pfm_ok ? 1 : 0, csv_ok ? 1 : 0,
                 json_ok ? 1 : 0);
    return 1;
  }

  if (base_config.frame_count > 1) {
    std::printf("Shadow frame %d/%d written at observer_time=%.6f:\n",
                frame_index + 1, base_config.frame_count, config.observer_time);
  } else {
    std::printf("Shadow outputs written to:\n");
  }
  std::printf("  %s\n", config.output_path);
  std::printf("  %s\n", config.hdr_output_path);
  std::printf("  %s\n", config.plot_output_path);
  std::printf("  %s\n", config.metadata_output_path);
  return 0;
}

int render_shadow_image() {
  const ShadowRenderConfig config = load_shadow_render_config();
  const FishboneMoncriefTorusConfig torus_config = load_fm_torus_config();
  if (std::abs(static_cast<double>(a)) > 1.0) {
    std::fprintf(stderr, "Shadow mode requires |spin| <= 1.\n");
    return 1;
  }
  if (!is_valid_render_config(config)) {
    std::fprintf(stderr, "Invalid shadow render configuration.\n");
    return 1;
  }
  if (!is_valid_fm_torus_config(torus_config)) {
    std::fprintf(stderr, "Invalid Fishbone-Moncrief torus configuration.\n");
    return 1;
  }

  FishboneMoncriefTorusModel torus_model;
  if (!build_fm_torus_model(static_cast<double>(a), torus_config,
                            torus_model)) {
    std::fprintf(stderr, "Failed to derive Fishbone-Moncrief torus model.\n");
    return 1;
  }

  for (int frame_index = 0; frame_index < config.frame_count; ++frame_index) {
    if (render_shadow_frame(config, torus_model, frame_index) != 0)
      return 1;
  }

  if (config.frame_count > 1) {
    const std::string input_pattern = build_sequence_pattern(config.output_path);
    const std::string video_output =
        build_sequence_video_output(config.output_path);
    std::printf("Sequence ready. Example ffmpeg command:\n");
    std::printf("  ffmpeg -framerate 24 -i %s -pix_fmt yuv420p %s\n",
                input_pattern.c_str(), video_output.c_str());
  }
  return 0;
}

} // namespace

int shadow_prob() { return render_shadow_image(); }
