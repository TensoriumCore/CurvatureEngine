#pragma once

#include <Geodesics.h>

#define DIM3 3
constexpr double DX = 0.09;
constexpr double DY = 0.09;
constexpr double DZ = 0.09;
#define NX 128
#define NY 128
#define NZ 128
#define GHOST 2  
#define NX_TOTAL (NX + 2*GHOST) 
#define NY_TOTAL (NY + 2*GHOST)
#define NZ_TOTAL (NZ + 2*GHOST)

using Matrix4x4 = std::array<std::array<double, NDIM>, NDIM>;
using Matrix3x3 = std::array<std::array<double, DIM3>, DIM3>;
using Vector3   = std::array<double, DIM3>;
using Vector4 = std::array<double, NDIM>;
using Tensor3D = std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>;
using Tensor4D  = std::array<std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>, DIM3>;
using Christoffel3D = std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>;
using Riemann3D = std::array<std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>, DIM3>;

class Grid {
    public:

		std::vector<double> dgammaX[3][3];
		std::vector<double> dgammaY[3][3];
		std::vector<double> dgammaZ[3][3];
		struct Cell2D {
			Matrix3x3 gamma;     
			Matrix3x3 gamma_inv; 
			Tensor3D Gamma3;    
			Matrix3x3 K;        
			Matrix3x3 Ricci;    
			double dgammaX[3][3];
			double dgammaY[3][3];
			double dgammaZ[3][3];
			double H;   
			double tilde_gamma[3][3];
			double chi;
			double momentum[3];
			double hamiltonian;
			double alpha;        
			Vector3 beta_cov;     
			Vector3 beta_con;  
			double tildeGamma[3];
			double Atilde[3][3];
			double beta[3];
			Matrix3x3 dgamma_dt;   
			Matrix3x3 dK_dt;      
			Vector3 dalpha_dx;    
			Tensor3D dK_dx;       
			Tensor3D dGamma3_dx;  
			double alphaStage[4];
			double betaStage[4][3];
			double Bstage[4][3];
			double rho;          
			Vector3 S_i;         
			Matrix3x3 S_ij;     
			double S; 
			double dgt[3][3]; 
			double dKt[3][3];  
			double gamma0[3][3];
			double K0[3][3];
			double alpha0;
			double beta0[3];
			double gammaStage[4][3][3]; 
			double Christoffel[3][3][3];
			double KStage[4][3][3];   
			double vx, vy, vz;
			double p;
			double T[4][4];
		};
		void compute_spectral_derivatives_for_gamma() ;
		void compute_energy_momentum_evolution(int i, int j, int k, double dt);
		void update_energy_momentum_tensor(int i, int j, int k);
		double compute_ricci_scalar(Grid &grid, int i, int j, int k);
		void export_fluid_slice(int j_slice);
		void initializeKerrData();
		void export_energy_momentum_tensor_slice(int slice_y);
		void update_fluid_velocity(int i, int j, int k, double dt);
		void compute_fluid_derivatives(int i, int j, int k);
		void initialize_grid();
		void evolve(Grid &grid_obj, double dtinitital, int nSteps);
		void extract_3p1(const Matrix4x4& g,
				const Matrix4x4& g_inv,  
				double* alpha,
				Vector3& beta_cov,
				Vector3& beta_con,
				Matrix3x3& gamma,
				Matrix3x3& gamma_inv);
		void initializeData(); 
		void compute_ricci_3d(
				Grid& grid_obj,  
				const Vector3& X,       
				const Tensor3D& Gamma3, 
				Matrix3x3& R3);
		void print_ricci_tensor(const Matrix3x3& R3);
		void calculate_christoffel_3D(const Vector3& X, Tensor3D& Gamma3, 
				const Matrix3x3& gamma, Matrix3x3 gamma_inv) ;
		void compute_ricci_3d(
				const Vector3& X,       
				const Tensor3D& Gamma3, 
				Matrix3x3& R3);
		double KUpAt(Grid &grid, int ip, int jp, int kp, int j_up, int i_low);
		void copyInitialState(Cell2D &cell);
		void updateIntermediateState(Cell2D &cell, double dtCoeff, int stageIndex);
		void storeStage(Cell2D &cell, int stage, double d_alpha_dt, double d_beta_dt[3]) ;
		Matrix3x3 compute_beta_gradient(int i, int j);
		void combineStages(Cell2D &cell, double dt);
		void initialize_grid(int Nr, int Ntheta, double r_min, double r_max, double theta_min, double theta_max);
		void export_vtk(const std::string& filename);
		double computeMaxSpeed();
		double computeCFL_dt(double CFL);
		void compute_constraints(Grid &grid_obj, int i, int j, int k, double &hamiltonian, double momentum[3]);
		void calculate_ricci_3d_from_riemann(const Riemann3D& Riemann, Matrix3x3& Ricci);
		void compute_time_derivatives(Grid &grid_obj, int i, int j, int k);
		void allocateGlobalGrid();
		void initializeData_Minkowski();
		void initializeKerrData(Grid &grid_obj);
		void compute_ricci_3D(Grid &grid_obj, int i, int j, int k, double Ricci[3][3]);
		double partialX_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);
		double partialY_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);
		double partialZ_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);

		double computeTraceK(Grid &grid, int i, int j, int k);
		double christoffelTerm(Grid &grid, int i, int j, int k, int i_comp);
		void compute_gauge_derivatives(Grid &grid_obj, int i, int j, int k, double &d_alpha_dt, double d_beta_dt[3]);
		
		Cell2D& getCell(int i, int j, int k) {
			return globalGrid[i][j][k];
		}

	private:
		std::vector<std::vector<std::vector<Grid::Cell2D>>> globalGrid;

};


double partialX_alpha(Grid &grid_obj, int i, int j, int k);
double partialY_alpha(Grid &grid_obj, int i, int j, int k);
double partialZ_alpha(Grid &grid_obj, int i, int j, int k);
double partialXX_alpha(Grid &grid_obj, int i, int j, int k);
double partialYY_alpha(Grid &grid_obj, int i, int j, int k);
double partialZZ_alpha(Grid &grid_obj, int i, int j, int k);
double partialX_betacomp(Grid &grid_obj, int i, int j, int k, int comp);
double partialY_betacomp(Grid &grid_obj, int i, int j, int k, int comp);
double partialZ_betacomp(Grid &grid_obj, int i, int j, int k, int comp);
double partialXY_alpha(Grid &grid_obj, int i, int j, int k);
double partialXZ_alpha(Grid &grid_obj, int i, int j, int k);
double partialYZ_alpha(Grid &grid_obj, int i, int j, int k);
double second_partial_alpha(Grid &grid_obj, int i, int j, int k, int a, int b);
bool invert_3x3(const double m[3][3], double inv[3][3]);
void apply_boundary_conditions(Grid &grid_obj);
void export_K_3D(Grid &grid_obj);
void export_alpha_slice(Grid &grid_obj, int j);
void export_gauge_slice(Grid &grid_obj, int j);
void export_K_slice(Grid &grid_obj, int j);
void export_gamma_slice(Grid &grid_obj, int j);
