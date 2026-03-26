#include "Connexion.h"
#include "Metric.h"
#include "Tensor.h"
#include "app/Problems.h"
#include "matrix.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

int Riemann_tensor(const char *metric) {
    Tensor tensor;
    Matrix matrix_obj;
    Connexion connexion;
    Metric metric_obj;
    
    float r0 = 6.0;
    std::array<float, NDIM> X = {0.0, r0, M_PI/2.0, 0.0};

    if (strcmp(metric, "kds") == 0) {
        printf("KDS metric calculation\n");
        metric_obj.calculate_metric_kds(X, metric_obj.gcov, metric_obj.gcon);
    } else if (strcmp(metric, "kerr-newman") == 0) {
        metric_obj.calculate_metric_kerr_newman(X, metric_obj.gcov, metric_obj.gcon);
    } else {
        metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);
        matrix_obj.print_matrix("metric_obj.gcov", metric_obj.gcov);
        matrix_obj.print_matrix("gcon", metric_obj.gcon);
    }

    connexion.calculate_christoffel(X, DELTA, connexion.Gamma, metric_obj.gcov, metric_obj.gcon, metric);

    for (int d = 0; d < NDIM; d++) {
        tensor.calculate_Gamma_at_offset(X, d, DELTA, DELTA, metric_obj.gcov, metric_obj.gcon, connexion.Gamma_plus_h[d], metric);
        tensor.calculate_Gamma_at_offset(X, d, -DELTA, DELTA, metric_obj.gcov, metric_obj.gcon, connexion.Gamma_minus_h[d], metric);
    }

    for (int d = 0; d < NDIM; d++) {
        tensor.calculate_Gamma_at_offset(X, d, DELTA/2.0, DELTA, metric_obj.gcov, metric_obj.gcon, connexion.Gamma_plus_half_h[d], metric);
        tensor.calculate_Gamma_at_offset(X, d, -DELTA/2.0, DELTA, metric_obj.gcov, metric_obj.gcon, connexion.Gamma_minus_half_h[d], metric);
    }

    tensor.calculate_riemann(connexion.Gamma, connexion.Gamma_plus_h,
                             connexion.Gamma_minus_h, connexion.Gamma_plus_half_h,
                             connexion.Gamma_minus_half_h, tensor.Riemann, DELTA);

    std::cout << "Riemann tensor compute finished" << std::endl;
    tensor.print_riemann(tensor.Riemann);
    tensor.contract_riemann(tensor.Riemann, tensor.Ricci, metric_obj.gcon);
    
    return 0;
}
