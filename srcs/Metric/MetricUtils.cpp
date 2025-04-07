#include <Geodesics.h>

void Metric::verify_metric(const std::array<std::array<float, NDIM>, NDIM>& g,
				const std::array<std::array<float, NDIM>, NDIM>& g_inv){
	
	Matrix matrix_obj;
    int i, j, k;
    float identity[NDIM][NDIM] = {0};
    float delta; 

    for (i = 0; i < NDIM; i++) {
        for (j = 0; j < NDIM; j++) {
            identity[i][j] = 0.0;
            for (k = 0; k < NDIM; k++) {
                identity[i][j] += gcon[i][k] * gcov[k][j];
            }
        }
    }

    for (i = 0; i < NDIM; i++) {
        for (j = 0; j < NDIM; j++) {
            if (i == j) {
                delta = 1.0;
            }
            else {
                delta = 0.0;
            }

            if (fabs(identity[i][j] - delta) > TOLERANCE) {
                printf("Erreur: identity[%d][%d] = %e with %e\n", i, j, identity[i][j], delta);
            }
        }
    }
	matrix_obj.check_inverse(gcov, gcon);
}
