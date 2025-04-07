#include <Geodesics.h>



void Tensor::calculate_Gamma_at_offset(const std::array<float, NDIM>& X, int direction, 
                                         float offset, float delta,
                                         Tensor::MatrixNDIM& gcov, 
                                         Tensor::MatrixNDIM& gcon, 
                                         Tensor::Christoffel3D& Gamma_slice, 
                                         const char* metric_type) {
    std::array<float, NDIM> X_offset = X;
    X_offset[direction] += offset;
    Tensor::Christoffel3D tempGamma{};
    Connexion connexion;
    Metric metric;
    
    if (strcmp(metric_type, "minkowski") == 0) {
        connexion.calculate_christoffel(X_offset, delta, tempGamma, gcov, gcon, "minkowski");
    } else if (strcmp(metric_type, "kerr") == 0 || strcmp(metric_type, "schwarzschild") == 0) {
        metric.calculate_metric(X_offset, gcov, gcon);
    }
    connexion.calculate_christoffel(X_offset, delta, tempGamma, gcov, gcon, metric_type);
    
    Gamma_slice = tempGamma;
}


float Tensor::richardson_derivative(
    const Tensor3D& Gamma_plus_h, 
    const Tensor3D& Gamma_minus_h,
    const Tensor3D& Gamma_plus_half_h,
    const Tensor3D& Gamma_minus_half_h,
    int rho, int mu, int nu, float h) 
{
    float diff_h = (Gamma_plus_h[rho][mu][nu] - Gamma_minus_h[rho][mu][nu]) / (2 * h);
    float diff_half_h = (Gamma_plus_half_h[rho][mu][nu] - Gamma_minus_half_h[rho][mu][nu]) / h;
    return (4 * diff_half_h - diff_h) / 3;
}

void Tensor::calculate_riemann(const Christoffel3D& Gamma, 
				const Christoffel4D& Gamma_plus_h, 
				const Christoffel4D& Gamma_minus_h, 
				const Christoffel4D& Gamma_plus_half_h, 
				const Christoffel4D& Gamma_minus_half_h,
				Riemann4D& Riemann, 
				float h) {
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    float dGamma_mu = richardson_derivative(
                        Gamma_plus_h[mu], Gamma_minus_h[mu], 
                        Gamma_plus_half_h[mu], Gamma_minus_half_h[mu], 
                        rho, nu, sigma, h);
                    float dGamma_nu = richardson_derivative(
                        Gamma_plus_h[nu], Gamma_minus_h[nu],
                        Gamma_plus_half_h[nu], Gamma_minus_half_h[nu],
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
	float Kretschmann_scalar = 0.0;
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    Kretschmann_scalar += Riemann[rho][sigma][mu][nu] * \
										  Riemann[rho][sigma][mu][nu];
                }
            }
        }
    }
    if (Kretschmann_scalar > 1e10) {
        printf("Kretschmann Scalar: INF (Singularity detected)\n");
    } else {
        printf("Kretschmann Scalar: %12.6f\n", Kretschmann_scalar);
    }
}
