#pragma once

#include <Geodesics.h>

#define DIM3 3
constexpr double DX = 0.08;
constexpr double DY = 0.08;
constexpr double DZ = 0.08;
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

struct alignas(32) Geometry {
    Matrix3x3 gamma;
    Matrix3x3 gamma_inv;
    Matrix3x3 Ricci;
    Matrix3x3 gamma0;
    double tilde_gamma[3][3];
    double tilde_gamma0[3][3];
    double tildgamma_inv[3][3];
    double dt_tilde_gamma[3][3];
    double dt_tildeGamma[3];
};

struct alignas(32) Connection {
	Tensor3D Gamma3;
	double tildeGamma[3];
	double Christoffel[3][3][3];
};

struct alignas(32) ExtrinsicCurvature {
    Matrix3x3 K;
    Matrix3x3 K0;
    double K_trace;
    double dt_K_trace;
    double H;
    double dKt[3][3];
};

struct alignas(32) AtildeVars {
    double Atilde[3][3];
    double Atilde0[3][3];
    double dt_Atilde[3][3];
    double AtildeStage[4][3][3];
};

struct alignas(32) Gauge {
    double alpha;
    double alpha0;
    double alphaStage[4];
    Vector3 dalpha_dx;

    double beta[3];
    double beta0[3];
    double betaStage[4][3];
    double Bstage[4][3];
	double dt_beta[3];
	double dt_alpha;
};

struct Matter {
    double rho;
    double momentum[3];
    double hamiltonian;
    double p;
    double vx, vy, vz;
    double T[4][4];
};


class Grid {
    public:

		std::vector<double> dgammaX[3][3];
		std::vector<double> dgammaY[3][3];
		std::vector<double> dgammaZ[3][3];
		double time = 0.0;
		struct alignas(32) Cell2D {
			Geometry geom;
			Connection conn;
			ExtrinsicCurvature curv;
			AtildeVars atilde;
			Gauge gauge;
			Matter matter;
			double t;
			double ADMmass;
			double chi;
			double dt_chi;
			double gammaStage[4][3][3];
			double KStage[4][3][3];
			double dgt[3][3];
		};

		struct alignas(32) GridStorage {
			Cell2D* cells;

			double* alpha;
			double* chi;
			double* Atilde[3][3];

			int nx, ny, nz;

			inline int idx(int i, int j, int k) const {
				return i * ny * nz + j * nz + k;
			}
		};

		void appendConstraintL2ToCSV(const std::string& filename, double time) const;
		void inject_BowenYork_Atilde(Grid &grid_obj, const Vector3 &P, const Vector3 &Coor);
		void logger_evolve(Grid &grid_obj, double dt, int nstep);
		double compute_ricci_scalar(Grid &grid, int i, int j, int k);
		void initialize_grid();
		void evolve(Grid &grid_obj, double dtinitital, int nSteps);
		double KUpAt(Grid &grid, int ip, int jp, int kp, int j_up, int i_low);
		void copyInitialState(Cell2D &cell);
		void updateIntermediateState(Cell2D &cell, double dtCoeff, int stageIndex);
		void storeStage(Cell2D &cell, int stage, double d_alpha_dt, double d_beta_dt[3]) ;
		void combineStages(Cell2D &cell, double dt);
		void initialize_grid(int Nr, int Ntheta, double r_min, double r_max, double theta_min, double theta_max);
		double computeMaxSpeed();
		double computeCFL_dt(double CFL);
		void compute_constraints(Grid &grid_obj, int i, int j, int k, double &hamiltonian, double momentum[3]);
		void compute_time_derivatives(Grid &grid_obj, int i, int j, int k);
		void allocateGlobalGrid();
		void initializeData_Minkowski();
		void initializeKerrData(Grid &grid_obj);
		void initializeBinaryKerrData(Grid &grid_obj);
		double partialX_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);
		double partialY_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);
		double partialZ_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);
		double computeTraceK(Grid &grid, int i, int j, int k);
		void compute_gauge_derivatives(Grid &grid_obj, int i, int j, int k, double &d_alpha_dt, double d_beta_dt[3]);
		void injectTTWave(Cell2D &cell, double x, double y, double z, double t);
		void solve_lichnerowicz(int max_iter, double tol, double dx, double dy, double dz);
		Cell2D& getCell(int i, int j, int k) {
			return globalGrid[i][j][k];
		}
		void export_Atildedt_slide(Grid &grid_obj, double time);
	private:
		std::vector<std::vector<std::vector<Grid::Cell2D>>> globalGrid;

};


double partialXX_alpha(Grid &grid_obj, int i, int j, int k);
double partialYY_alpha(Grid &grid_obj, int i, int j, int k);
double partialZZ_alpha(Grid &grid_obj, int i, int j, int k);
double second_partial_alpha(Grid &grid_obj, int i, int j, int k, int a, int b);
bool invert_3x3(const double m[3][3], double inv[3][3]);
void apply_boundary_conditions(Grid &grid_obj);
void export_K_3D(Grid &grid_obj);
void export_alpha_slice(Grid &grid_obj, int j);
void export_gauge_slice(Grid &grid_obj, int j);
void export_K_slice(Grid &grid_obj, int j);
void export_gamma_slice(Grid &grid_obj, int j, double time);
void export_tilde_gamma_3D(Grid &grid_obj);
