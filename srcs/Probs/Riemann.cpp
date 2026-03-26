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
    using Christoffel3D = Connexion::Christoffel3D;
    using Christoffel4D = Connexion::Christoffel4D;
    Tensor tensor;
    Matrix matrix_obj;
    Connexion connexion;
    Metric metric_obj;
    constexpr float riemann_delta = 5.0e-2f;
    Christoffel3D gamma{};
    Christoffel4D gamma_plus_h{};
    Christoffel4D gamma_minus_h{};
    Tensor::Riemann4D riemann{};
    Tensor::MatrixNDIM ricci{};
    
    float r0 = 6.0;
    std::array<float, NDIM> X = {0.0, r0, M_PI/2.0, 0.0};

    if (strcmp(metric, "kds") == 0) {
        printf("KDS metric calculation\n");
    }
    metric_obj.calculate_metric_by_name(metric, X, metric_obj.gcov, metric_obj.gcon);
    matrix_obj.print_matrix("metric_obj.gcov", metric_obj.gcov);
    matrix_obj.print_matrix("gcon", metric_obj.gcon);

    connexion.calculate_christoffel(X, riemann_delta, gamma, metric_obj.gcov, metric_obj.gcon, metric);

    for (int d = 0; d < NDIM; d++) {
        tensor.calculate_Gamma_at_offset(X, d, riemann_delta, riemann_delta, metric_obj.gcov, metric_obj.gcon, gamma_plus_h[d], metric);
        tensor.calculate_Gamma_at_offset(X, d, -riemann_delta, riemann_delta, metric_obj.gcov, metric_obj.gcon, gamma_minus_h[d], metric);
    }

    tensor.calculate_riemann(gamma, gamma_plus_h, gamma_minus_h, riemann, riemann_delta);
    const double kretschmann =
        tensor.calculate_kretschmann_scalar(riemann, metric_obj.gcov,
                                            metric_obj.gcon);
    tensor.contract_riemann(riemann, ricci);
    const double ricci_scalar = tensor.calculate_ricci_scalar(ricci, metric_obj.gcon);

    printf("Kretschmann Scalar: %.12f\n", kretschmann);
    std::cout << "Riemann tensor compute finished" << std::endl;
    tensor.print_riemann(riemann);
    matrix_obj.print_matrix("Ricci", ricci);
    printf("Ricci Scalar: %12.6f\n", ricci_scalar);
    
    return 0;
}
