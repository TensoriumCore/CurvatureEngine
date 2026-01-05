#include <Geodesics.h>

extern float (*geodesic_points)[5];
extern int num_points;
extern float a;

using namespace curvatureengine::simd; // Indispensable pour load/store/broadcast

namespace {

struct LightChristoffelContext {
  Connexion *connexion;
  Metric *metric;
  const char *metric_name;
  float step;
};

void evaluate_light_christoffel(const double coords[NDIM],
                                ChristoffelTensor gamma, void *ctx) {
  auto *context = static_cast<LightChristoffelContext *>(ctx);
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

void evaluate_light_christoffel_adapter(const VEC_TYPE coords[NDIM],
                                        ChristoffelTensorVec& gamma, void *ctx) {
    double pos_unpacked[NDIM][4];
    for(int i = 0; i < NDIM; ++i) {
        store(pos_unpacked[i], coords[i]);
    }

    double gamma_buffer[NDIM][NDIM][NDIM][4];

    for(int lane = 0; lane < 4; ++lane) {
        double local_coords[NDIM];
        for(int i = 0; i < NDIM; ++i) {
            local_coords[i] = pos_unpacked[i][lane];
        }

        ChristoffelTensor local_gamma;
        evaluate_light_christoffel(local_coords, local_gamma, ctx);

        for(int i = 0; i < NDIM; ++i)
            for(int j = 0; j < NDIM; ++j)
                for(int k = 0; k < NDIM; ++k)
                    gamma_buffer[i][j][k][lane] = local_gamma[i][j][k];
    }

    for(int i = 0; i < NDIM; ++i) {
        for(int j = 0; j < NDIM; ++j) {
            for(int k = 0; k < NDIM; ++k) {
                gamma[i][j][k] = load(gamma_buffer[i][j][k]);
            }
        }
    }
}

} // namespace

int light_geodesics_prob() {
  Connexion connexion;
  Metric metric_obj;
  float r0 = 100.0;
  std::array<float, NDIM> X = {0.0, r0, M_PI / 4.0, 0.0};
  
  metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);
  float g_tt = metric_obj.gcov[0][0];
  float g_tphi = metric_obj.gcov[0][3];
  float g_phiphi = metric_obj.gcov[3][3];
  float Omega = 1.0 / (pow(r0, 1.5) + a);
  float denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
  float vt = 1.0 / sqrt(fabs(denom));
  float v[NDIM] = {vt, 0.0, 0.0, 3.5f * Omega * vt};
  float norm =
      g_tt * v[0] * v[0] + 2.0 * g_tphi * v[0] * v[3] + g_phiphi * v[3] * v[3];
  
  float dt = 0.00910;
  LightChristoffelContext ctx{&connexion, &metric_obj, "kerr", DELTA};
  
  VEC_TYPE X_avx[NDIM], v_avx[NDIM];
  for (int i = 0; i < NDIM; i++) {
    X_avx[i] = broadcast(X[i]);
    v_avx[i] = broadcast(v[i]);
  }

  auto start = std::chrono::high_resolution_clock::now();
  
  geodesic_AVX(X_avx, v_avx, max_dt + 4, evaluate_light_christoffel_adapter, &ctx,
               broadcast(dt));
               
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<float> elapsed_seconds = end - start;
  printf("Elapsed time: %f\n", elapsed_seconds.count());
  write_vtk_file("output/light_geodesic.vtk");
  if (geodesic_points != NULL) {
    free(geodesic_points);
  }
  return 0;
}
