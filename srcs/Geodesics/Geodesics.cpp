#include "app/RuntimeState.h"
#include "core/Constants.h"
#include "core/GeodesicIntegrator.h"

#include <algorithm>
#include <cmath>

using namespace curvatureengine::simd;

namespace {

struct StageDerivativesVec {
  VEC_TYPE dx[NDIM];
  VEC_TYPE dv[NDIM];
};

void compute_vector_stage(const VEC_TYPE state_x[NDIM],
                          const VEC_TYPE state_v[NDIM],
                          ChristoffelEvalFnVec eval_fn, void *ctx,
                          StageDerivativesVec &stage) {
  ChristoffelTensorVec gamma;
  if (eval_fn) {
    eval_fn(state_x, gamma, ctx);
  }

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

void integrate_geodesic_vector(VEC_TYPE x[NDIM], VEC_TYPE v[NDIM],
                               float lambda_max, ChristoffelEvalFnVec eval_fn,
                               void *ctx, VEC_TYPE step_size) {

  VEC_TYPE lambda = zero();
  double l_max = static_cast<double>(lambda_max);
  double r_horizon_limit = 1.0 + std::sqrt(1.0 - a * a) + 0.15;

  while (lane0(lambda) < l_max) {
    if (lane0(x[1]) < r_horizon_limit) {
      break;
    }

    double h_val = std::min(lane0(step_size), l_max - lane0(lambda));
    VEC_TYPE h = broadcast(h_val);

    StageDerivativesVec k1, k2, k3, k4;
    VEC_TYPE tmp_x[NDIM], tmp_v[NDIM];

    compute_vector_stage(x, v, eval_fn, ctx, k1);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k1.dx[i] * h * 0.5;
      tmp_v[i] = v[i] + k1.dv[i] * h * 0.5;
    }

    compute_vector_stage(tmp_x, tmp_v, eval_fn, ctx, k2);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k2.dx[i] * h * 0.5;
      tmp_v[i] = v[i] + k2.dv[i] * h * 0.5;
    }

    compute_vector_stage(tmp_x, tmp_v, eval_fn, ctx, k3);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k3.dx[i] * h;
      tmp_v[i] = v[i] + k3.dv[i] * h;
    }

    compute_vector_stage(tmp_x, tmp_v, eval_fn, ctx, k4);
    VEC_TYPE h6 = h * (1.0 / 6.0);
    for (int i = 0; i < NDIM; ++i) {
      x[i] =
          x[i] + h6 * (k1.dx[i] + 2.0 * k2.dx[i] + 2.0 * k3.dx[i] + k4.dx[i]);
      v[i] =
          v[i] + h6 * (k1.dv[i] + 2.0 * k2.dv[i] + 2.0 * k3.dv[i] + k4.dv[i]);
    }

    lambda = lambda + h;
    store_geodesic_point_AVX(x, static_cast<float>(lane0(lambda)));
  }
}

VEC_TYPE integrate_raytrace(VEC_TYPE x[NDIM], VEC_TYPE v[NDIM],
                            float lambda_max, ChristoffelEvalFnVec eval_fn,
                            void *ctx, VEC_TYPE step_size) {
  VEC_TYPE lambda = zero();
  double l_max = static_cast<double>(lambda_max);
  double r_horizon_limit = 1.0 + std::sqrt(1.0 - a * a) + 0.15;
  double r_escape_limit = 250.0;
  alignas(32) double r_check[4];

  while (lane0(lambda) < l_max) {
    store(r_check, x[1]);
    bool active = false;
    for (int i = 0; i < 4; ++i) {
      if (r_check[i] > r_horizon_limit && r_check[i] < r_escape_limit) {
        active = true;
        break;
      }
    }
    if (!active) {
      break;
    }

    double h_val = std::min(lane0(step_size), l_max - lane0(lambda));
    VEC_TYPE h = broadcast(h_val);

    StageDerivativesVec k1, k2, k3, k4;
    VEC_TYPE tmp_x[NDIM], tmp_v[NDIM];

    compute_vector_stage(x, v, eval_fn, ctx, k1);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k1.dx[i] * h * 0.5;
      tmp_v[i] = v[i] + k1.dv[i] * h * 0.5;
    }
    compute_vector_stage(tmp_x, tmp_v, eval_fn, ctx, k2);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k2.dx[i] * h * 0.5;
      tmp_v[i] = v[i] + k2.dv[i] * h * 0.5;
    }
    compute_vector_stage(tmp_x, tmp_v, eval_fn, ctx, k3);
    for (int i = 0; i < NDIM; ++i) {
      tmp_x[i] = x[i] + k3.dx[i] * h;
      tmp_v[i] = v[i] + k3.dv[i] * h;
    }
    compute_vector_stage(tmp_x, tmp_v, eval_fn, ctx, k4);

    VEC_TYPE h6 = h * (1.0 / 6.0);
    for (int i = 0; i < NDIM; ++i) {
      x[i] =
          x[i] + h6 * (k1.dx[i] + 2.0 * k2.dx[i] + 2.0 * k3.dx[i] + k4.dx[i]);
      v[i] =
          v[i] + h6 * (k1.dv[i] + 2.0 * k2.dv[i] + 2.0 * k3.dv[i] + k4.dv[i]);
    }
    lambda = lambda + h;
  }
  return x[1];
}
} // namespace

void geodesic_AVX(VEC_TYPE x[NDIM], VEC_TYPE v[NDIM], float lambda_max,
                  ChristoffelEvalFnVec evaluator, void *ctx,
                  VEC_TYPE step_size) {
  integrate_geodesic_vector(x, v, lambda_max, evaluator, ctx, step_size);
}

VEC_TYPE geodesic_raytrace_AVX(VEC_TYPE x[NDIM], VEC_TYPE v[NDIM],
                               float lambda_max, ChristoffelEvalFnVec evaluator,
                               void *ctx, VEC_TYPE step_size) {
  return integrate_raytrace(x, v, lambda_max, evaluator, ctx, step_size);
}
