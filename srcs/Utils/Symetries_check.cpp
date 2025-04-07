#include <Geodesics.h>

void Tensor::check_riemann_symmetries(const Riemann4D& Riemann, float tolerance) {
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    if (fabs(Riemann[rho][sigma][mu][nu] + Riemann[rho][sigma][nu][mu]) > tolerance) {
                        printf("Non-antisymmetry R[%d][%d][%d][%d] and R[%d][%d][%d][%d]\n", rho, sigma, mu, nu, rho, sigma, nu, mu);
                    }
                    if (fabs(Riemann[rho][sigma][mu][nu] - Riemann[sigma][rho][mu][nu]) > tolerance) {
                        printf("Non-symmetry R[%d][%d][%d][%d] and R[%d][%d][%d][%d]\n", rho, sigma, mu, nu, sigma, rho, mu, nu);
                    }
                    if (fabs(Riemann[rho][sigma][mu][nu] + Riemann[rho][mu][nu][sigma] + Riemann[rho][nu][sigma][mu]) > tolerance) {
                        printf("Bianchi identity violation at R[%d][%d][%d][%d]\n", rho, sigma, mu, nu);
                    }
                }
            }
        }
    }
}
void Connexion::check_symmetry_christoffel(const Christoffel3D& gamma) {
    printf("\nChristoffel symbols verif with delta = %e:\n", TOLERANCE);
    for (int lambda = 0; lambda < NDIM; lambda++) {
        for (int mu = 0; mu < NDIM; mu++) {
            for (int nu = mu; nu < NDIM; nu++) { 
                float diff = fabs(gamma[lambda][mu][nu] - gamma[lambda][nu][mu]);
                if (diff > TOLERANCE) {
                    printf("!symmetry at Gamma^%d_%d%d: |%f - %f| = %f\n", 
                           lambda, mu, nu, gamma[lambda][mu][nu], gamma[lambda][nu][mu], diff);
                }
            }
        }
    }
    printf("Symmetric !\n");
}

