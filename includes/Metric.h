#pragma once	

#include "Geodesics.h"

class Metric {
	public:
		std::array<std::array<float, NDIM>, NDIM> gcov{};
		std::array<std::array<float, NDIM>, NDIM> gcon{};
		std::array<std::array<float, NDIM>, NDIM> gcovK{};
		std::array<std::array<float, NDIM>, NDIM> gconK{};
		std::array<std::array<float, NDIM>, NDIM> gconMinkowski{};
		std::array<std::array<float, NDIM>, NDIM> gcov_half{};
		std::array<std::array<float, NDIM>, NDIM> gcon_half{};
		std::array<std::array<float, NDIM>, NDIM> gcovMinkowski{};

		Metric() = default;

		~Metric() = default;


		
		
		void calculate_metric(const std::array<float, NDIM>& x, 
				std::array<std::array<float, NDIM>, NDIM>& g,
				std::array<std::array<float, NDIM>, NDIM>& g_inv);
		void verify_metric(const std::array<std::array<float, NDIM>, NDIM>& g,
				const std::array<std::array<float, NDIM>, NDIM>& g_inv);
		void calculate_metric_kds(const std::array<float, NDIM>& x, 
				std::array<std::array<float, NDIM>, NDIM>& g,
				std::array<std::array<float, NDIM>, NDIM>& g_inv);
		void calculate_metric_kerr_newman(const std::array<float, NDIM>& x, 
				std::array<std::array<float, NDIM>, NDIM>& g,
				std::array<std::array<float, NDIM>, NDIM>& g_inv);

};
