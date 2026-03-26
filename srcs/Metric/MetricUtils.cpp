#include "Metric.h"
#include "matrix.h"

#include <cmath>
#include <cstring>
#include <cstdio>

void Metric::calculate_metric_by_name(
    const char* metric_name,
    const std::array<float, NDIM>& x,
    std::array<std::array<float, NDIM>, NDIM>& g,
    std::array<std::array<float, NDIM>, NDIM>& g_inv) {
    if (metric_name != nullptr && std::strcmp(metric_name, "kds") == 0) {
        calculate_metric_kds(x, g, g_inv);
        return;
    }
    if (metric_name != nullptr && std::strcmp(metric_name, "kerr-newman") == 0) {
        calculate_metric_kerr_newman(x, g, g_inv);
        return;
    }

    calculate_metric(x, g, g_inv);
}

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
                identity[i][j] += g_inv[i][k] * g[k][j];
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
	matrix_obj.check_inverse(g, g_inv);
}
