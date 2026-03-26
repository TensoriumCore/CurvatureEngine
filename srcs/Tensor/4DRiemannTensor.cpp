#include "Connexion.h"
#include "Metric.h"
#include "Tensor.h"

#include <cstring>

void Tensor::calculate_Gamma_at_offset(const std::array<float, NDIM>& X, int direction, 
                                         float offset, float delta,
                                         const Tensor::MatrixNDIM& gcov, 
                                         const Tensor::MatrixNDIM& gcon, 
                                         Tensor::Christoffel3D& Gamma_slice, 
                                         const char* metric_type) {
    std::array<float, NDIM> X_offset = X;
    X_offset[direction] += offset;
    Tensor::Christoffel3D tempGamma{};
    Connexion connexion;
    Metric metric;
    Tensor::MatrixNDIM gcov_local = gcov;
    Tensor::MatrixNDIM gcon_local = gcon;
    
    if (strcmp(metric_type, "minkowski") != 0) {
        metric.calculate_metric_by_name(metric_type, X_offset, gcov_local, gcon_local);
    }
    connexion.calculate_christoffel(X_offset, delta, tempGamma, gcov_local, gcon_local, metric_type);
    
    Gamma_slice = tempGamma;
}

float Tensor::central_difference_derivative(
    const Tensor3D& Gamma_plus_h, 
    const Tensor3D& Gamma_minus_h,
    int rho, int mu, int nu, float h) 
{
    float diff_h = (Gamma_plus_h[rho][mu][nu] - Gamma_minus_h[rho][mu][nu]) / (2 * h);
    return diff_h;
}

void Tensor::calculate_riemann(const Christoffel3D& Gamma, 
				const Christoffel4D& Gamma_plus_h, 
				const Christoffel4D& Gamma_minus_h, 
				Riemann4D& Riemann, 
				float h) {
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    float dGamma_mu = central_difference_derivative(
                        Gamma_plus_h[mu], Gamma_minus_h[mu], 
                        rho, nu, sigma, h);
                    float dGamma_nu = central_difference_derivative(
                        Gamma_plus_h[nu], Gamma_minus_h[nu],
                        rho, mu, sigma, h);

                    float Gamma_terms = 0.0;
                    for (int lambda = 0; lambda < NDIM; lambda++) {
                        Gamma_terms += Gamma[rho][mu][lambda] * Gamma[lambda][nu][sigma]
                                     - Gamma[rho][nu][lambda] * Gamma[lambda][mu][sigma];
                    }
                    Riemann[rho][sigma][mu][nu] = dGamma_mu - dGamma_nu + Gamma_terms;
                }
            }
        }
    }
}

double Tensor::calculate_kretschmann_scalar(const Riemann4D& Riemann,
                                            const MatrixNDIM& g_cov,
                                            const MatrixNDIM& g_inv) {
    double kretschmann = 0.0;

    for (int rho = 0; rho < NDIM; ++rho) {
        for (int sigma = 0; sigma < NDIM; ++sigma) {
            for (int mu = 0; mu < NDIM; ++mu) {
                for (int nu = 0; nu < NDIM; ++nu) {
                    const double lhs = static_cast<double>(Riemann[rho][sigma][mu][nu]);
                    for (int alpha = 0; alpha < NDIM; ++alpha) {
                        for (int beta = 0; beta < NDIM; ++beta) {
                            for (int gamma = 0; gamma < NDIM; ++gamma) {
                                for (int delta = 0; delta < NDIM; ++delta) {
                                    kretschmann +=
                                        static_cast<double>(g_cov[rho][alpha]) *
                                        static_cast<double>(g_inv[sigma][beta]) *
                                        static_cast<double>(g_inv[mu][gamma]) *
                                        static_cast<double>(g_inv[nu][delta]) *
                                        lhs *
                                        static_cast<double>(Riemann[alpha][beta][gamma][delta]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return kretschmann;
}
