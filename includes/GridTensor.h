#pragma once 

#include <Geodesics.h>

class GridTensor {
	public:
		GridTensor() = default;
		~GridTensor() = default;    
		friend class Grid;
		void export_christoffel_slice(Grid &grid_obj, int j);
		void compute_extrinsic_curvature(Grid &grid_obj, int i, int j, int k, \
											 double dx, double dy, double dz);
		void compute_Atilde(Grid &grid_obj, int i, int j, int k);
	protected:
		void compute_christoffel_3D(Grid &grid_obj, int i, int j, int k, double christof[3][3][3]);
		void compute_dt_tildeGamma(Grid &grid_obj, int i, int j, int k, double dt_tildeGamma[3]); 
		void compute_tildeGamma(Grid &grid_obj, int i, int j, int k, double tildeGamma[3]);
		void compute_partial_christoffel(Grid &grid_obj, int i, int j, int k, int dim, double partialGamma[3][3][3][3], double d);
		void compute_ricci_conformal_factor(Grid &grid_obj, int i, int j, int k, double RicciChi[3][3]);
		void compute_ricci_BSSN(Grid &grid_obj, int i, int j, int k, double Ricci[3][3]);
		void compute_ricci_3D_conformal(Grid &grid_obj, int i, int j, int k, double Ricci[3][3]);
		void compute_ricci_3d(
				Grid& grid_obj,  
				const Vector3& X,       
				const Tensor3D& Gamma3, 
				Matrix3x3& R3);
		double partialX_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialY_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialZ_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialX_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialY_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialZ_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
		double partialX_gammaSpec(Grid &grid_obj, int i, int j, int k, int a, int b);
		double compute_momentum_i(Grid &grid, int i, int j, int k, int i_comp);
		double christoffelTerm(Grid &grid, int i, int j, int k, int i_comp);
};


class BSSNevolve {
	public:
		BSSNevolve() = default;
		~BSSNevolve() = default;
		friend class Grid;
		void compute_dt_chi(Grid &grid_obj, int i, int j, int k, double &dt_chi);
		void compute_dt_tilde_gamma(Grid &grid_obj, int i, int j, int k, double dt_tg[3][3]);

	protected:

};
