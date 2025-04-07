#include <Geodesics.h>
extern float (*geodesic_points)[5];
extern int num_points;
extern float a;

int light_geodesics_prob() {
	Connexion connexion;
	Metric metric_obj;
	float r0 = 100.0;
    std::array<float, NDIM> X = {0.0, r0, M_PI/4.0, 0.0};;
	metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);
	float g_tt = metric_obj.gcov[0][0];
	float g_tphi = metric_obj.gcov[0][3];
	float g_phiphi = metric_obj.gcov[3][3];
	float Omega = 1.0 / (pow(r0, 1.5) + a);
	float denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
	float vt = 1.0 / sqrt(fabs(denom));
	float v[NDIM] = {vt, 0.0, 0.0, 3.5f * Omega * vt};
	float norm = g_tt * v[0] * v[0] + 2.0 * g_tphi * v[0] * v[3] + g_phiphi * v[3] * v[3];
	float dt = 0.00910;
	float christoffel[NDIM][NDIM][NDIM];
	connexion.calculate_christoffel(X, DELTA, connexion.Gamma, metric_obj.gcov, metric_obj.gcon, "kerr");
	__m256d X_avx[NDIM], v_avx[NDIM];
	for (int i = 0; i < NDIM; i++) {
		X_avx[i] = _mm256_set1_pd(X[i]);
		v_avx[i] = _mm256_set1_pd(v[i]);
	}
	__m256d christoffel_avx[NDIM][NDIM][NDIM];
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			for (int k = 0; k < NDIM; k++) {
				christoffel_avx[i][j][k] = _mm256_set1_pd(connexion.Gamma[i][j][k]);
			}
		}
	}

	auto start = std::chrono::high_resolution_clock::now();
	geodesic_AVX(X_avx, v_avx, max_dt + 4, ( __m256d (*)[NDIM][NDIM] )christoffel_avx, _mm256_set1_pd(dt));
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<float> elapsed_seconds = end - start;
	printf("Elapsed time: %f\n", elapsed_seconds.count());
	write_vtk_file("output/light_geodesic.vtk");
	if (geodesic_points != NULL) {
		free(geodesic_points);
	}
	return 0;
}
