#pragma once

#include "core/Types.h"

class Tensor {
    public:
        using Christoffel3D =
            std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>;
        using Christoffel4D = std::array<
            std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>, NDIM>;
        using Riemann4D = std::array<
            std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>, NDIM>;
        using MatrixNDIM = std::array<std::array<float, NDIM>, NDIM>;
        using Tensor3D =
            std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>;

        Tensor() = default;

        void calculate_Gamma_at_offset(const std::array<float, NDIM>& X,
                                       int direction, float offset, float delta,
                                       const MatrixNDIM& gcov,
                                       const MatrixNDIM& gcon,
                                       Christoffel3D& Gamma_slice,
                                       const char* metric_type);
        float central_difference_derivative(const Tensor3D& Gamma_plus_h,
                                            const Tensor3D& Gamma_minus_h,
                                            int rho, int mu, int nu, float h);
        void calculate_riemann(const Christoffel3D& Gamma,
                               const Christoffel4D& Gamma_plus_h,
                               const Christoffel4D& Gamma_minus_h,
                               Riemann4D& Riemann, float h);

        void print_riemann(const Riemann4D& Riemann);
        void check_riemann_symmetries(const Riemann4D& Riemann, float tolerance);
        void contract_riemann(const Riemann4D& Riemann, MatrixNDIM& Ricci);
        double calculate_ricci_scalar(const MatrixNDIM& Ricci,
                                      const MatrixNDIM& g_inv);
        double calculate_kretschmann_scalar(const Riemann4D& Riemann,
                                            const MatrixNDIM& g_cov,
                                            const MatrixNDIM& g_inv);
};
