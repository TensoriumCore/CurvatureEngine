#pragma once 

#include <Geodesics.h>

class GridTensor {
	public:
		GridTensor() = default;
		~GridTensor() = default;    
		friend class Grid;
		void export_christoffel_slice(Grid &grid_obj, int j);
		void compute_extrinsic_curvature(Grid &grid_obj, int i, int j, int k, \
											 float dx, float dy, float dz);
		void compute_Atilde(Grid &grid_obj, int i, int j, int k);
	protected:
		void compute_christoffel_3D(Grid &grid_obj, int i, int j, int k, float christof[3][3][3]);
		void compute_dt_tildeGamma(Grid &grid_obj, int i, int j, int k, float dt_tildeGamma[3]); 
		void compute_tildeGamma(Grid &grid_obj, int i, int j, int k, float tildeGamma[3]);
		void compute_partial_christoffel(Grid &grid_obj, int i, int j, int k, int dim, float partialGamma[3][3][3][3], float d);
		void compute_ricci_conformal_factor(Grid &grid_obj, int i, int j, int k, float RicciChi[3][3]);
		void compute_ricci_BSSN(Grid &grid_obj, int i, int j, int k, float Ricci[3][3]);
		void compute_ricci_3D_conformal(Grid &grid_obj, int i, int j, int k, float Ricci[3][3]);
		void compute_ricci_3d(
				Grid& grid_obj,  
				const Vector3& X,       
				const Tensor3D& Gamma3, 
				Matrix3x3& R3);
		float partialX_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		float partialY_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		float partialZ_gamma(Grid &grid_obj, int i, int j, int k, int a, int b);
		float partialX_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
		float partialY_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
		float partialZ_Kij(Grid &grid_obj, int i, int j, int k, int a, int b);
		float partialX_gammaSpec(Grid &grid_obj, int i, int j, int k, int a, int b);
		float compute_momentum_i(Grid &grid, int i, int j, int k, int i_comp);
		float christoffelTerm(Grid &grid, int i, int j, int k, int i_comp);
};


class BSSNevolve {
	public:
		BSSNevolve() = default;
		~BSSNevolve() = default;
		friend class Grid;
		void compute_dt_chi(Grid &grid_obj, int i, int j, int k, float &dt_chi);
		void compute_dt_tilde_gamma(Grid &grid_obj, int i, int j, int k, float dt_tg[3][3]);

	protected:

};
