#include "Metric.h"
#include "app/RuntimeState.h"
#include "matrix.h"

#include <cmath>

void Metric::calculate_metric_kds(const std::array<float, NDIM>& x, 
                                  std::array<std::array<float, NDIM>, NDIM>& g,
                                  std::array<std::array<float, NDIM>, NDIM>& g_inv) {
    Matrix matrix_obj;
    float r = x[1];
    float theta = x[2];
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);

    float Sigma = r * r + a * a * cos_theta * cos_theta;
    float Delta_r = (1.0 - (KERR_DE_SITTER_LAMBDA * r * r) / 3.0f) * (r * r + a * a) - 2.0f * M * r;
    float Delta_theta = 1.0f + (KERR_DE_SITTER_LAMBDA * a * a / 3.0f) * cos_theta * cos_theta;
    float Xi = 1.0f + (KERR_DE_SITTER_LAMBDA * a * a) / 3.0f;
    float sin_theta2 = sin_theta * sin_theta;
    float r2_plus_a2 = r * r + a * a;

    for (auto& row : g) {
        row.fill(0.0);
    }
    for (auto& row : g_inv) {
        row.fill(0.0);
    }

    g[0][0] = (-Delta_r + a * a * Delta_theta * sin_theta2) / Sigma;
    g[1][1] = Sigma / Delta_r;
    g[2][2] = Sigma / Delta_theta;
    g[3][3] = (sin_theta2 / (Sigma * Xi * Xi)) *
              (r2_plus_a2 * r2_plus_a2 * Delta_theta -
               a * a * sin_theta2 * Delta_r);
    g[0][3] = (a * sin_theta2 / (Sigma * Xi)) *
              (Delta_r - r2_plus_a2 * Delta_theta);
    g[3][0] = g[0][3];

    matrix_obj.inverse_matrix(g, g_inv);

}
