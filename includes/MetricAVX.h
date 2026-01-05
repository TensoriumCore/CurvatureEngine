#pragma once

#include "Geodesics.h"
#include "SimdConfig.h"
#include <cmath>

extern float a;

namespace MetricSIMD {
using namespace curvatureengine::simd;

inline VEC_TYPE sin_vec(VEC_TYPE x) {
  alignas(32) double buf[4];
  store(buf, x);
  for (int i = 0; i < 4; ++i)
    buf[i] = std::sin(buf[i]);
  return load(buf);
}

inline VEC_TYPE cos_vec(VEC_TYPE x) {
  alignas(32) double buf[4];
  store(buf, x);
  for (int i = 0; i < 4; ++i)
    buf[i] = std::cos(buf[i]);
  return load(buf);
}

inline void calculate_metric_avx(const VEC_TYPE x[NDIM], VEC_TYPE g[NDIM][NDIM],
                                 VEC_TYPE g_inv[NDIM][NDIM]) {

  VEC_TYPE one = broadcast(1.0);
  VEC_TYPE two = broadcast(2.0);
  VEC_TYPE m_vec = broadcast(M);
  VEC_TYPE a_vec = broadcast(a);

  VEC_TYPE r = x[1];
  VEC_TYPE theta = x[2];
  VEC_TYPE st = sin_vec(theta);
  VEC_TYPE ct = cos_vec(theta);
  VEC_TYPE st2 = st * st;
  VEC_TYPE ct2 = ct * ct;
  VEC_TYPE r2 = r * r;
  VEC_TYPE a2 = a_vec * a_vec;

  VEC_TYPE Sigma = r2 + a2 * ct2;
  VEC_TYPE Delta = r2 - (two * m_vec * r) + a2;

  VEC_TYPE epsilon = broadcast(1e-6);
  VEC_TYPE Delta_safe = Delta + epsilon;
  VEC_TYPE Sigma_safe = Sigma + epsilon;

  VEC_TYPE inv_Sigma = one / Sigma_safe;
  VEC_TYPE inv_Delta = one / Delta_safe;
  VEC_TYPE z = zero();
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      g[i][j] = z;
      g_inv[i][j] = z;
    }

  g[0][0] = -(one - (two * m_vec * r * inv_Sigma));
  g[1][1] = Sigma * inv_Delta;
  g[2][2] = Sigma;
  VEC_TYPE term_mix = (two * m_vec * r * a_vec * st2) * inv_Sigma;
  g[0][3] = -term_mix;
  g[3][0] = g[0][3];
  g[3][3] = st2 * (r2 + a2 + (a_vec * term_mix));

  g_inv[1][1] = Delta * inv_Sigma;
  g_inv[2][2] = inv_Sigma;
  VEC_TYPE det_block = g[0][0] * g[3][3] - g[0][3] * g[0][3];
  VEC_TYPE inv_det = one / det_block;

  g_inv[0][0] = g[3][3] * inv_det;
  g_inv[3][3] = g[0][0] * inv_det;
  g_inv[0][3] = -g[0][3] * inv_det;
  g_inv[3][0] = g_inv[0][3];
}

} // namespace MetricSIMD

namespace {
using namespace curvatureengine::simd;

inline void evaluate_christoffel_native_avx(const VEC_TYPE x[NDIM],
                                            ChristoffelTensorVec &gamma,
                                            void *ctx) {
  const float step_val = 1e-4f;
  const VEC_TYPE h = broadcast(step_val);
  const VEC_TYPE inv_2h = broadcast(0.5f / step_val);
  const VEC_TYPE half = broadcast(0.5);
  const VEC_TYPE zero_v = zero();

  VEC_TYPE g_up[NDIM][NDIM];
  VEC_TYPE g_unused[NDIM][NDIM];

  MetricSIMD::calculate_metric_avx(x, g_unused, g_up);

  VEC_TYPE dg[NDIM][NDIM][NDIM];
  VEC_TYPE x_shifted[NDIM];
  for (int i = 0; i < NDIM; i++)
    x_shifted[i] = x[i];
  for (int sigma = 0; sigma < NDIM; ++sigma) {
    VEC_TYPE original_x = x_shifted[sigma];
    x_shifted[sigma] = original_x + h;
    VEC_TYPE g_fwd[NDIM][NDIM], g_inv_trash[NDIM][NDIM];
    MetricSIMD::calculate_metric_avx(x_shifted, g_fwd, g_inv_trash);

    x_shifted[sigma] = original_x - h;
    VEC_TYPE g_bwd[NDIM][NDIM];
    MetricSIMD::calculate_metric_avx(x_shifted, g_bwd, g_inv_trash);

    x_shifted[sigma] = original_x;

    for (int mu = 0; mu < NDIM; ++mu) {
      for (int nu = 0; nu < NDIM; ++nu) {
        dg[mu][nu][sigma] = (g_fwd[mu][nu] - g_bwd[mu][nu]) * inv_2h;
      }
    }
  }

  for (int k = 0; k < NDIM; ++k) {
    for (int i = 0; i < NDIM; ++i) {
      for (int j = 0; j < NDIM; ++j) {

        VEC_TYPE sum = zero_v;
        for (int l = 0; l < NDIM; ++l) {
          VEC_TYPE term = dg[l][i][j] + dg[l][j][i] - dg[i][j][l];
          sum = fmadd(g_up[k][l], term, sum);
        }

        gamma[k][i][j] = sum * half;
      }
    }
  }
}
} // namespace
