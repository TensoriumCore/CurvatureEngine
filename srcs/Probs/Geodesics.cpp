#include <Geodesics.h>
extern float (*geodesic_points)[5];
extern int num_points;
extern float a;

namespace {

struct ChristoffelEvalContext {
  Connexion *connexion;
  Metric *metric;
  const char *metric_name;
  float step;
};

void evaluate_christoffel_at_point(const double coords[NDIM],
                                   ChristoffelTensor gamma, void *ctx) {
  auto *context = static_cast<ChristoffelEvalContext *>(ctx);
  Connexion::VectorNDIM Xf{};
  for (int i = 0; i < NDIM; ++i) {
    Xf[i] = static_cast<float>(coords[i]);
  }
  Connexion::Christoffel3D gamma_float{};
  Connexion::MatrixNDIM g_local{};
  Connexion::MatrixNDIM g_inv_local{};
  context->metric->calculate_metric(Xf, g_local, g_inv_local);
  context->connexion->calculate_christoffel(Xf, context->step, gamma_float,
                                            g_local, g_inv_local,
                                            context->metric_name);
  for (int mu = 0; mu < NDIM; ++mu) {
    for (int alpha = 0; alpha < NDIM; ++alpha) {
      for (int beta = 0; beta < NDIM; ++beta) {
        gamma[mu][alpha][beta] =
            static_cast<double>(gamma_float[mu][alpha][beta]);
      }
    }
  }
}

} // namespace

int Geodesics_prob() {
  Connexion connexion;
  Metric metric_obj;
  float r0 = 19.3;
  std::array<float, NDIM> X = {0.0, r0, M_PI / 3.8, 0.0};
  ;
  metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);
  float g_tt = metric_obj.gcov[0][0];
  float g_tphi = metric_obj.gcov[0][3];
  float g_phiphi = metric_obj.gcov[3][3];
  float Omega = 1.0 / (pow(r0, 1.5) + a);
  float denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
  float vt = 1.0 / sqrt(denom);
  float v[NDIM] = {vt, 0.0, 0.0, Omega * vt};
  float norm =
      g_tt * v[0] * v[0] + 2.0 * g_tphi * v[0] * v[3] + g_phiphi * v[3] * v[3];
  printf("Norme quadrivecteur = %e (doit être proche de 0)\n", norm);
  float dt = 0.0910;
  ChristoffelEvalContext ctx{&connexion, &metric_obj, "kerr", DELTA};

  VEC_TYPE X_avx[NDIM], v_avx[NDIM];
  for (int i = 0; i < NDIM; i++) {
    X_avx[i] = curvatureengine::simd::broadcast(X[i]);
    v_avx[i] = curvatureengine::simd::broadcast(v[i]);
  }

  printf("begin geodesic\n");
  auto start = std::chrono::high_resolution_clock::now();
  geodesic_AVX(X_avx, v_avx, max_dt + 4, evaluate_christoffel_at_point, &ctx,
    curvatureengine::simd::broadcast(dt));
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<float> elapsed_seconds = end - start;
  printf("Elapsed time: %f\n", elapsed_seconds.count());
  write_vtk_file("Output/geodesic.vtk");
  if (geodesic_points != NULL) {
    free(geodesic_points);
  }
  return 0;
}
