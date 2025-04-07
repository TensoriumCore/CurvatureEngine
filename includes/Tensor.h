#pragma once

#include "Geodesics.h"

class Tensor {
	public:
		using Christoffel3D = std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>;
		using Christoffel4D = std::array<std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>, NDIM>;
		using Riemann4D = std::array<std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>, NDIM>;
		using MatrixNDIM = std::array<std::array<float, NDIM>, NDIM>;
		using VectorNDIM = std::array<float, NDIM>;
		using Tensor3D = std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>;
		
		Christoffel3D gamma{};
		Christoffel3D Gamma_plus_h{};
		Christoffel3D Gamma_minus_h{};
		Christoffel4D Gamma_plus_half_h{};
		Christoffel4D Gamma_minus_half_h{};
        Riemann4D Riemann{};  
        MatrixNDIM Ricci{};
		MatrixNDIM Ricci_scalar{};
		MatrixNDIM Einstein{};
		MatrixNDIM Einstein_tensor{};
			
		Tensor() = default;
	
		void calculate_Gamma_at_offset(const std::array<float, NDIM>& X, int direction, 
				float offset, float delta,
				MatrixNDIM& gcov, 
				MatrixNDIM& gcon, 
				Christoffel3D& Gamma_slice, 
				const char* metric_type);
		float richardson_derivative(
				const Tensor3D& Gamma_plus_h, 
				const Tensor3D& Gamma_minus_h,
				const Tensor3D& Gamma_plus_half_h,
				const Tensor3D& Gamma_minus_half_h,
				int rho, int mu, int nu, float h); 
		void calculate_riemann(const Christoffel3D& Gamma, 
				const Christoffel4D& Gamma_plus_h, 
				const Christoffel4D& Gamma_minus_h, 
				const Christoffel4D& Gamma_plus_half_h, 
				const Christoffel4D& Gamma_minus_half_h,
				Riemann4D& Riemann, 
				float h);

		void initialize_riemann_tensor(Riemann4D& R);
		void print_riemann(const Riemann4D& Riemann);
		void check_riemann_symmetries(const Riemann4D& Riemann, float tolerance);
		void contract_riemann(const Riemann4D& Riemann, MatrixNDIM& Ricci, const MatrixNDIM& g_inv);
};

