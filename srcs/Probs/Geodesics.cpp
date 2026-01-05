#include "MetricAVX.h"
#include <Geodesics.h>

extern float (*geodesic_points)[5];
extern int num_points;
extern float a;

using namespace curvatureengine::simd;

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

void evaluate_christoffel_adapter(const VEC_TYPE coords[NDIM],
                                  ChristoffelTensorVec &gamma, void *ctx) {
  double pos_unpacked[NDIM][4];

  for (int i = 0; i < NDIM; ++i) {
    store(pos_unpacked[i], coords[i]);
  }

  double gamma_buffer[NDIM][NDIM][NDIM][4];

  for (int lane = 0; lane < 4; ++lane) {
    double local_coords[NDIM];
    for (int i = 0; i < NDIM; ++i) {
      local_coords[i] = pos_unpacked[i][lane];
    }

    ChristoffelTensor local_gamma;
    evaluate_christoffel_at_point(local_coords, local_gamma, ctx);

    for (int i = 0; i < NDIM; ++i)
      for (int j = 0; j < NDIM; ++j)
        for (int k = 0; k < NDIM; ++k)
          gamma_buffer[i][j][k][lane] = local_gamma[i][j][k];
  }

  for (int i = 0; i < NDIM; ++i) {
    for (int j = 0; j < NDIM; ++j) {
      for (int k = 0; k < NDIM; ++k) {
        gamma[i][j][k] = load(gamma_buffer[i][j][k]);
      }
    }
  }
}

} // namespace

extern std::atomic<int> global_idx;

using namespace curvatureengine::simd;

int Geodesics_prob() {
  Connexion connexion;
  Metric metric_obj;

  const size_t max_points = 10000000;
  if (geodesic_points)
    free(geodesic_points);

  if (posix_memalign((void **)&geodesic_points, 32,
                     max_points * sizeof(*geodesic_points)) != 0) {
    fprintf(stderr, "Allocation failed\n");
    return -1;
  }

  global_idx.store(0);

  float r0 = 19.3;
  std::array<float, NDIM> X = {0.0, r0, M_PI / 3.8, 0.0};

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
  printf("Norme quadrivecteur = %e\n", norm);

  float dt = 0.0910;

  VEC_TYPE X_avx[NDIM], v_avx[NDIM];
  for (int i = 0; i < NDIM; i++) {
    X_avx[i] = broadcast(X[i]);
    v_avx[i] = broadcast(v[i]);
  }

  printf("begin geodesic\n");
  auto start = std::chrono::high_resolution_clock::now();

  geodesic_AVX(X_avx, v_avx, max_dt + 4, evaluate_christoffel_native_avx,
               nullptr, broadcast(dt));

  auto end = std::chrono::high_resolution_clock::now();

  num_points = global_idx.load();

  std::chrono::duration<float> elapsed_seconds = end - start;
  printf("Elapsed time: %f\n", elapsed_seconds.count());
  printf("Points computed: %d\n", num_points);

  write_vtk_file("Output/geodesic.vtk");

  if (geodesic_points != NULL) {
    free(geodesic_points);
    geodesic_points = NULL;
  }

  return 0;
}
