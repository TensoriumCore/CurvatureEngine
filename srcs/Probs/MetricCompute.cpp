#include "Metric.h"
#include "app/Problems.h"
#include "core/Constants.h"
#include "matrix.h"

#include <array>
#include <cmath>
#include <cstdio>

int Metric_prob() {
	Metric metric_obj;
	Matrix matrix_obj;
	float r0 = 20.0;
    std::array<float, NDIM> X = {0.0, r0, M_PI/4.0, 0.0};;
	printf("=====================================================\n");
	metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);
	printf("Kerr metric calculated\n");
	matrix_obj.print_matrix("g", metric_obj.gcov);
	matrix_obj.print_matrix("g_inv", metric_obj.gcon);

	printf("=====================================================\n");
	std::array<std::array<float, NDIM>, NDIM> gcov_KN, gcon_KN;	
	metric_obj.calculate_metric_kerr_newman(X, gcov_KN, gcon_KN);
	printf("Kerr-Newman metric calculated\n");
	matrix_obj.print_matrix("g", gcov_KN);
	matrix_obj.print_matrix("g_inv", gcon_KN);

	printf("=====================================================\n");
	std::array<std::array<float, NDIM>, NDIM> gcov_kds, gcon_kds;
	metric_obj.calculate_metric_kds(X, gcov_kds, gcon_kds);
	printf("Kerr-de Sitter metric calculated\n");
	matrix_obj.print_matrix("g", gcov_kds);
	matrix_obj.print_matrix("g_inv", gcon_kds);
	return 0;
}
