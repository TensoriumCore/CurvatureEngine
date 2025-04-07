#include "Connexion.h"
#include <Geodesics.h>

void initialize_riemann_tensor(float R[NDIM][NDIM][NDIM][NDIM]) {
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    R[mu][nu][rho][sigma] = 0.0; 
                }
            }
        }
    }
}


void Tensor::print_riemann(const Riemann4D& Riemann) {
    const float threshold = 1e-10; 
    printf("\nRiemann Tensor (Non-zero components):\n"); 
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    float value = Riemann[rho][sigma][mu][nu];
                    if (fabs(value) > threshold) {
                        printf("Riemann[%d][%d][%d][%d] = %12.6f\n", rho, sigma, mu, nu, value);
                    }
                }
            }
        }
    }
}

void Connexion::print_christoffel_matrix(const Christoffel3D& gamma) {
    printf("\nChristoffel Symbols:\n");
    for (int lambda = 0; lambda < NDIM; lambda++) {
        printf("\nGamma^%d:\n", lambda);
        for (int mu = 0; mu < NDIM; mu++) {
            for (int nu = 0; nu < NDIM; nu++) {
                printf("%12.6f\t", gamma[lambda][mu][nu]);
            }
            printf("\n");
        }
    }
}
