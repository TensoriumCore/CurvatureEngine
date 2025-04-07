#include <Geodesics.h>

extern float (*geodesic_points)[5];
extern int num_points;
extern float a;
float Lambda = 1e-4;

void Metric::calculate_metric_kds(const std::array<float, NDIM>& x, 
                                  std::array<std::array<float, NDIM>, NDIM>& g,
                                  std::array<std::array<float, NDIM>, NDIM>& g_inv) {
    Matrix matrix_obj;
    Metric metric_obj;
    float r = x[1];
    float theta = x[2];
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);

    float Sigma = r * r + a * a * cos_theta * cos_theta;
    float Delta_r = (1.0 - (Lambda * r * r) / 3.0) * (r * r + a * a) - 2.0 * M * r;
    float Delta_theta = 1.0 + (Lambda * a * a / 3.0) * cos_theta * cos_theta;
    float Xi = 1.0 - (Lambda * a * a) / 3.0;

    for (auto& row : g) {
        row.fill(0.0);
    }
    for (auto& row : g_inv) {
        row.fill(0.0);
    }

    g[0][0] = - (Delta_r / (Sigma * Xi * Xi));
    g[1][1] = Sigma / Delta_r;
    g[2][2] = Sigma / Delta_theta;
    g[3][3] = (sin_theta * sin_theta / (Sigma * Xi * Xi)) * (r * r + a * a) * (r * r + a * a);
    g[0][3] = - (2.0 * M * r * a * sin_theta * sin_theta) / (Sigma * Xi * Xi);
    g[3][0] = g[0][3];

    matrix_obj.inverse_matrix(g, g_inv);

    if (a == 0.0 && Lambda == 0.0)
        printf("Schwarzschild metric calculated\n");
    else if (a != 0.0 && Lambda == 0.0)
        printf("Kerr metric calculated\n");
    else 
        printf("Kerr-de Sitter metric calculated\n");

    matrix_obj.print_matrix("g", g);
}


