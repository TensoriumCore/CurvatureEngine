
#include <Geodesics.h>

extern float (*geodesic_points)[5];
extern int num_points;
extern float a;

void Metric::calculate_metric(const std::array<float, NDIM>& x, 
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
    float Delta = r * r - 2.0 * M * r + a * a;
    
    g[0][0] = -(1.0 - (2.0 * M * r) / Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = sin_theta2 * (r * r + a * a + (2.0 * M * r * a * a * sin_theta2) / Sigma);
    g[0][3] = - (2.0 * M * r * a * sin_theta2) / Sigma;
    g[3][0] = g[0][3];

    matrix_obj.inverse_matrix(g, g_inv);
}


