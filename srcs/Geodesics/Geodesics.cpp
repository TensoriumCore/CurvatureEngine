#include <Geodesics.h>
#include <algorithm>

extern float a;

namespace {

struct StageDerivatives {
  double dx[NDIM];
  double dv[NDIM];
};

void compute_scalar_stage(const double state_x[NDIM],
                          const double state_v[NDIM], ChristoffelEvalFn eval_fn,
                          void *ctx, StageDerivatives &stage) {
  ChristoffelTensor gamma{};
  if (eval_fn) {
    eval_fn(state_x, gamma, ctx);
  }

  for (int mu = 0; mu < NDIM; ++mu) {
    stage.dx[mu] = state_v[mu];
    double acc = 0.0;
    if (eval_fn) {
      for (int alpha = 0; alpha < NDIM; ++alpha) {
        for (int beta = 0; beta < NDIM; ++beta) {
          acc -= gamma[mu][alpha][beta] * state_v[alpha] * state_v[beta];
        }
      }
    }
    stage.dv[mu] = acc;
  }
}

void store_geodesic_point_scalar(const double x_scalar[NDIM], double lambda) {
  using curvatureengine::simd::broadcast;

  VEC_TYPE x_vec[NDIM];
  for (int mu = 0; mu < NDIM; ++mu) {
    x_vec[mu] = broadcast(x_scalar[mu]);
  }
  store_geodesic_point_AVX(x_vec, static_cast<float>(lambda));
}

void integrate_geodesic_scalar(double x_scalar[NDIM], double v_scalar[NDIM],
                               double lambda_max, ChristoffelEvalFn eval_fn,
                               void *ctx, double step_size) {
  double lambda = 0.0;
  while (lambda < lambda_max) {
    const double h =
        std::min(step_size, static_cast<double>(lambda_max) - lambda);

    StageDerivatives k1{}, k2{}, k3{}, k4{};
    double temp_x[NDIM];
    double temp_v[NDIM];

    compute_scalar_stage(x_scalar, v_scalar, eval_fn, ctx, k1);

    for (int mu = 0; mu < NDIM; ++mu) {
      temp_x[mu] = x_scalar[mu] + 0.5 * h * k1.dx[mu];
      temp_v[mu] = v_scalar[mu] + 0.5 * h * k1.dv[mu];
    }
    compute_scalar_stage(temp_x, temp_v, eval_fn, ctx, k2);

    for (int mu = 0; mu < NDIM; ++mu) {
      temp_x[mu] = x_scalar[mu] + 0.5 * h * k2.dx[mu];
      temp_v[mu] = v_scalar[mu] + 0.5 * h * k2.dv[mu];
    }
    compute_scalar_stage(temp_x, temp_v, eval_fn, ctx, k3);

    for (int mu = 0; mu < NDIM; ++mu) {
      temp_x[mu] = x_scalar[mu] + h * k3.dx[mu];
      temp_v[mu] = v_scalar[mu] + h * k3.dv[mu];
    }
    compute_scalar_stage(temp_x, temp_v, eval_fn, ctx, k4);

    for (int mu = 0; mu < NDIM; ++mu) {
      const double sum_dx =
          k1.dx[mu] + 2.0 * k2.dx[mu] + 2.0 * k3.dx[mu] + k4.dx[mu];
      const double sum_dv =
          k1.dv[mu] + 2.0 * k2.dv[mu] + 2.0 * k3.dv[mu] + k4.dv[mu];
      x_scalar[mu] += (h / 6.0) * sum_dx;
      v_scalar[mu] += (h / 6.0) * sum_dv;
    }

    lambda += h;
    store_geodesic_point_scalar(x_scalar, lambda);
  }
}

} // namespace

void geodesic_AVX(VEC_TYPE x[4], VEC_TYPE v[4], float lambda_max,
                  ChristoffelEvalFn evaluator, void *ctx, VEC_TYPE step_size) {
  using curvatureengine::simd::broadcast;
  using curvatureengine::simd::lane0;

  double x_scalar[NDIM];
  double v_scalar[NDIM];
  for (int mu = 0; mu < NDIM; ++mu) {
    x_scalar[mu] = lane0(x[mu]);
    v_scalar[mu] = lane0(v[mu]);
  }

  const double h = lane0(step_size);
  integrate_geodesic_scalar(x_scalar, v_scalar, lambda_max, evaluator, ctx, h);

  for (int mu = 0; mu < NDIM; ++mu) {
    x[mu] = broadcast(x_scalar[mu]);
    v[mu] = broadcast(v_scalar[mu]);
  }
}
