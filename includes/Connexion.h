#pragma once 

#include "core/Types.h"

class Connexion {
	public:
		using Christoffel3D = std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>;
		using Christoffel4D = std::array<std::array<std::array<std::array<float, NDIM>, NDIM>, NDIM>, NDIM>;
		using MatrixNDIM = std::array<std::array<float, NDIM>, NDIM>;
		using VectorNDIM = std::array<float, NDIM>;

		void calculate_christoffel(const VectorNDIM& X, float h,
				Christoffel3D& gamma,
				std::array<std::array<float, NDIM>, NDIM>& g,
				std::array<std::array<float, NDIM>, NDIM>& g_inv, 
				const char* metric); 

		void check_symmetry_christoffel(const Christoffel3D& gamma);
		void print_christoffel(const Christoffel3D& Gamma);
		void print_christoffel_matrix(const Christoffel3D& gamma);
};
