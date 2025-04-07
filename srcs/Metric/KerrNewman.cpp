
#include <Geodesics.h>

extern float (*geodesic_points)[5];
extern int num_points;
extern float a;
float Q = 0.9;

void Metric::calculate_metric_kerr_newman(const std::array<float, NDIM>& x, 
				std::array<std::array<float, NDIM>, NDIM>& g,
				std::array<std::array<float, NDIM>, NDIM>& g_inv) {
	Matrix matrix_obj;
	Metric metric_obj;
	float r = x[1];
    float theta = x[2];
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    float sin_theta2 = sin_theta * sin_theta;
    float cos_theta2 = cos_theta * cos_theta;

    float Sigma = r * r + a * a * cos_theta2;
    float Delta = r * r - 2.0 * M * r + a * a + Q * Q;

	for (auto& row : g) {
		row.fill(0.0);
	} 
	for (auto& row : g_inv) {
		row.fill(0.0);
	}

    g[0][0] = -(1.0 - (2.0 * M * r - Q * Q) / Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = sin_theta2 * (r * r + a * a + (2.0 * M * r - Q * Q) * a * a * sin_theta2 / Sigma);
    g[0][3] = -((2.0 * M * r - Q * Q) * a * sin_theta2) / Sigma;
    g[3][0] = g[0][3];

    matrix_obj.inverse_matrix(g, g_inv);

    if (a == 0.0 && Q == 0.0) 
        printf("Schwarzschild metric calculated\n");
    else if (Q == 0.0) 
        printf("Kerr metric calculated\n");
	else
        printf("Kerr-Newman metric calculated\n");

    matrix_obj.print_matrix("g", g);
    matrix_obj.print_matrix("g_inv", g_inv);
}



