#include "Tensor.h"

void Tensor::contract_riemann(const Riemann4D& Riemann, MatrixNDIM& Ricci) {
    for (auto &row : Ricci) {
        row.fill(0.0);
    }
    
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                Ricci[mu][nu] += Riemann[rho][mu][rho][nu];
            }
        }
    }
}

double Tensor::calculate_ricci_scalar(const MatrixNDIM& Ricci, const MatrixNDIM& g_inv) {
    double Ricci_scalar = 0.0;
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            Ricci_scalar += g_inv[mu][nu] * Ricci[mu][nu];
        }
    }
    return Ricci_scalar;
}
