#include "Connexion.h"
#include "Metric.h"

#include <cmath>
#include <cstring>
#include <iostream>

namespace {

using MatrixNDIMDouble = std::array<std::array<double, NDIM>, NDIM>;
using Christoffel3DDouble =
    std::array<std::array<std::array<double, NDIM>, NDIM>, NDIM>;

void compute_schwarzschild_metric(const std::array<double, NDIM> &coords,
                                  Connexion::MatrixNDIM &gout,
                                  Connexion::MatrixNDIM &ginv_out,
                                  MatrixNDIMDouble &gout_d,
                                  MatrixNDIMDouble &ginv_d) {
  const double r = coords[1];
  const double theta = coords[2];
  const double sin_theta = std::sin(theta);
  const double sin_theta2 = sin_theta * sin_theta;
  const double f = 1.0 - (2.0 / r);
  const double inv_f = 1.0 / f;

  gout = {};
  ginv_out = {};
  gout[0][0] = static_cast<float>(-f);
  gout[1][1] = static_cast<float>(inv_f);
  gout[2][2] = static_cast<float>(r * r);
  gout[3][3] = static_cast<float>(r * r * sin_theta2);

  ginv_out[0][0] = static_cast<float>(-1.0 / f);
  ginv_out[1][1] = static_cast<float>(f);
  ginv_out[2][2] = static_cast<float>(1.0 / (r * r));
  ginv_out[3][3] = static_cast<float>(1.0 / (r * r * sin_theta2));

  gout_d = {};
  ginv_d = {};
  gout_d[0][0] = -f;
  gout_d[1][1] = inv_f;
  gout_d[2][2] = r * r;
  gout_d[3][3] = r * r * sin_theta2;

  ginv_d[0][0] = -1.0 / f;
  ginv_d[1][1] = f;
  ginv_d[2][2] = 1.0 / (r * r);
  ginv_d[3][3] = 1.0 / (r * r * sin_theta2);
}

void populate_metric(const char *metric_name, Metric &metric_obj,
                     const Connexion::VectorNDIM &coords,
                     Connexion::MatrixNDIM &gout,
                     Connexion::MatrixNDIM &ginv_out) {
  if (std::strcmp(metric_name, "kds") == 0) {
    metric_obj.calculate_metric_kds(coords, gout, ginv_out);
    return;
  }
  if (std::strcmp(metric_name, "kerr-newman") == 0) {
    metric_obj.calculate_metric_kerr_newman(coords, gout, ginv_out);
    return;
  }

  metric_obj.calculate_metric(coords, gout, ginv_out);
}

} // namespace

void Connexion::calculate_christoffel(
    const VectorNDIM &X, float h, Christoffel3D &gamma,
    std::array<std::array<float, NDIM>, NDIM> &g,
    std::array<std::array<float, NDIM>, NDIM> &g_inv, const char *metric) {
  Christoffel3DDouble dg{};
  Christoffel3DDouble tmp_double{};
  VectorNDIM Xh = X, Xl = X;
  MatrixNDIM gh{}, gl{};
  MatrixNDIM g_inv_h{}, g_inv_l{};
  MatrixNDIM g_base{};
  MatrixNDIM g_inv_base{};
  MatrixNDIMDouble g_base_d{};
  MatrixNDIMDouble g_inv_base_d{};
  MatrixNDIMDouble gh_d{}, gl_d{};
  MatrixNDIMDouble ginv_tmp_h{}, ginv_tmp_l{};
  std::array<double, NDIM> X_double{};
  for (int i = 0; i < NDIM; ++i) {
    X_double[i] = static_cast<double>(X[i]);
  }
  Metric metric_obj;

  const float step = (h == 0.0f) ? DELTA : h;
  const double denom = 2.0 * static_cast<double>(step);

  auto compute_metric = [&](const VectorNDIM &coords, MatrixNDIM &gout,
                            MatrixNDIM &ginv_out, MatrixNDIMDouble &gout_d,
                            MatrixNDIMDouble &ginv_d) {
    if (strcmp(metric, "minkowski") == 0) {
      gout = g_base;
      ginv_out = g_inv_base;
      gout_d = g_base_d;
      ginv_d = g_inv_base_d;
    } else {
      populate_metric(metric, metric_obj, coords, gout, ginv_out);
      for (int i = 0; i < NDIM; ++i) {
        for (int j = 0; j < NDIM; ++j) {
          gout_d[i][j] = static_cast<double>(gout[i][j]);
          ginv_d[i][j] = static_cast<double>(ginv_out[i][j]);
        }
      }
    }
  };

  if (strcmp(metric, "minkowski") == 0) {
    g_base = g;
    g_inv_base = g_inv;
    for (int i = 0; i < NDIM; ++i) {
      for (int j = 0; j < NDIM; ++j) {
        g_base_d[i][j] = static_cast<double>(g_base[i][j]);
        g_inv_base_d[i][j] = static_cast<double>(g_inv_base[i][j]);
      }
    }
  } else if (strcmp(metric, "schwarzschild_test") == 0) {
    compute_schwarzschild_metric(X_double, g_base, g_inv_base, g_base_d,
                                 g_inv_base_d);
  } else {
    populate_metric(metric, metric_obj, X, g_base, g_inv_base);
    for (int i = 0; i < NDIM; ++i) {
      for (int j = 0; j < NDIM; ++j) {
        g_base_d[i][j] = static_cast<double>(g_base[i][j]);
        g_inv_base_d[i][j] = static_cast<double>(g_inv_base[i][j]);
      }
    }
  }

  g = g_base;
  g_inv = g_inv_base;

  for (int mu = 0; mu < NDIM; mu++) {
    if (strcmp(metric, "schwarzschild_test") == 0) {
      auto Xh_d = X_double;
      auto Xl_d = X_double;
      Xh_d[mu] += step;
      Xl_d[mu] -= step;
      compute_schwarzschild_metric(Xh_d, gh, g_inv_h, gh_d, ginv_tmp_h);
      compute_schwarzschild_metric(Xl_d, gl, g_inv_l, gl_d, ginv_tmp_l);
    } else {
      Xh = X;
      Xl = X;
      Xh[mu] += step;
      Xl[mu] -= step;

      if (strcmp(metric, "minkowski") == 0) {
        gh = g_base;
        gl = g_base;
        gh_d = g_base_d;
        gl_d = g_base_d;
      } else {
        compute_metric(Xh, gh, g_inv_h, gh_d, ginv_tmp_h);
        compute_metric(Xl, gl, g_inv_l, gl_d, ginv_tmp_l);
      }
    }

    for (int lam = 0; lam < NDIM; lam++) {
      for (int nu = 0; nu < NDIM; nu++) {
        double diff = gh_d[lam][nu] - gl_d[lam][nu];
        dg[lam][nu][mu] = diff / denom;
      }
    }
  }

  for (int lam = 0; lam < NDIM; lam++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int mu = 0; mu < NDIM; mu++) {
        tmp_double[lam][nu][mu] =
            0.5 * (dg[nu][lam][mu] + dg[mu][lam][nu] - dg[mu][nu][lam]);
      }
    }
  }

  for (int lam = 0; lam < NDIM; lam++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int mu = 0; mu < NDIM; mu++) {
        double sum = 0.0;
        for (int kap = 0; kap < NDIM; kap++) {
          sum += g_inv_base_d[lam][kap] * tmp_double[kap][nu][mu];
        }
        gamma[lam][nu][mu] = static_cast<float>(sum);
      }
    }
  }

#ifdef CURVATUREENGINE_DEBUG_CHRISTOFFEL
  std::cout << "Christoffel symbols calculated\n";
  print_christoffel_matrix(gamma);
  check_symmetry_christoffel(gamma);
#endif
}
