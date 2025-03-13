#pragma once 

#include <Geodesics.h>


template <class T>
inline double fourth_order_diff(const T &plus2, const T &plus1, 
                         const T &minus1, const T &minus2, double dx)
{
    return (-plus2 + 8.0*plus1 - 8.0*minus1 + minus2) / (12.0*dx);
}


template <class T>
inline double second_order_diff(const T &plus1, const T &minus1, double dx)
{
    return (plus1 - minus1)/(2.0*dx);
}

class GridTensor {
	public:
		GridTensor() = default;
		~GridTensor() = default;    
		friend class Grid;
		void export_christoffel_slice(Grid &grid_obj, int j);
		void compute_extrinsic_curvature(Grid &grid_obj, int i, int j, int k, \
											 double dx, double dy, double dz);
	protected:
		void compute_christoffel_3D(Grid &grid_obj, int i, int j, int k, double christof[3][3][3]);
		void compute_dt_tildeGamma(int i, int j, int k, double dt_tildeGamma[3]); 
		void compute_tildeGamma(int i, int j, int k, double tildeGamma[3]);
		void compute_partial_christoffel(Grid &grid_obj, int i, int j, int k, int dim, double partialGamma[3][3][3][3], double d);
		double partialX_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialY_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialZ_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialX_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialY_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialZ_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
};

