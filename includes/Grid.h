#pragma once

#include <Geodesics.h>

#define DIM3 3
constexpr float DX = 0.08;
constexpr float DY = 0.08;
constexpr float DZ = 0.08;
#define NX 128
#define NY 128
#define NZ 128
#define GHOST 2  
#define NX_TOTAL (NX + 2*GHOST) 
#define NY_TOTAL (NY + 2*GHOST)
#define NZ_TOTAL (NZ + 2*GHOST)

using Matrix4x4 = std::array<std::array<float, NDIM>, NDIM>;
using Matrix3x3 = std::array<std::array<float, DIM3>, DIM3>;
using Vector3   = std::array<float, DIM3>;
using Vector4 = std::array<float, NDIM>;
using Tensor3D = std::array<std::array<std::array<float, DIM3>, DIM3>, DIM3>;
using Tensor4D  = std::array<std::array<std::array<std::array<float, DIM3>, DIM3>, DIM3>, DIM3>;
using Christoffel3D = std::array<std::array<std::array<float, DIM3>, DIM3>, DIM3>;
using Riemann3D = std::array<std::array<std::array<std::array<float, DIM3>, DIM3>, DIM3>, DIM3>;

struct alignas(32) Geometry {
    Matrix3x3 gamma;
    Matrix3x3 gamma_inv;
    Matrix3x3 Ricci;
    Matrix3x3 gamma0;
    float tilde_gamma[3][3];
    float tilde_gamma0[3][3];
    float tildgamma_inv[3][3];
    float dt_tilde_gamma[3][3];
    float dt_tildeGamma[3];
};

struct alignas(32) Connection {
	Tensor3D Gamma3;
	float tildeGamma[3];
	float Christoffel[3][3][3];
};

struct alignas(32) ExtrinsicCurvature {
    Matrix3x3 K;
    Matrix3x3 K0;
    float K_trace;
    float dt_K_trace;
    float H;
    float dKt[3][3];
};

struct alignas(32) AtildeVars {
    float Atilde[3][3];
    float Atilde0[3][3];
    float dt_Atilde[3][3];
    float AtildeStage[4][3][3];
};

struct alignas(32) Gauge {
    float alpha;
    float alpha0;
    float alphaStage[4];
    Vector3 dalpha_dx;

    float beta[3];
    float beta0[3];
    float betaStage[4][3];
    float Bstage[4][3];
	float dt_beta[3];
	float dt_alpha;
};

struct Matter {
    float rho;
    float momentum[3];
    float hamiltonian;
    float p;
    float vx, vy, vz;
    float T[4][4];
};


class Grid {
    public:

		std::vector<float> dgammaX[3][3];
		std::vector<float> dgammaY[3][3];
		std::vector<float> dgammaZ[3][3];
		float time = 0.0;
		struct alignas(32) Cell2D {
			Geometry geom;
			Connection conn;
			ExtrinsicCurvature curv;
			AtildeVars atilde;
			Gauge gauge;
			Matter matter;
			float t;
			float ADMmass;
			float chi;
			float dt_chi;
			float gammaStage[4][3][3];
			float KStage[4][3][3];
			float dgt[3][3];
		};

		struct alignas(32) GridStorage {
			Cell2D* cells;

			float* alpha;
			float* chi;
			float* Atilde[3][3];

			int nx, ny, nz;

			inline int idx(int i, int j, int k) const {
				return i * ny * nz + j * nz + k;
			}
		};

		void appendConstraintL2ToCSV(const std::string& filename, float time) const;
		void inject_BowenYork_Atilde(Grid &grid_obj, const Vector3 &P, const Vector3 &Coor);
		void logger_evolve(Grid &grid_obj, float dt, int nstep);
		float compute_ricci_scalar(Grid &grid, int i, int j, int k);
		void initialize_grid();
		void evolve(Grid &grid_obj, float dtinitital, int nSteps);
		float KUpAt(Grid &grid, int ip, int jp, int kp, int j_up, int i_low);
		void copyInitialState(Cell2D &cell);
		void updateIntermediateState(Cell2D &cell, float dtCoeff, int stageIndex);
		void storeStage(Cell2D &cell, int stage, float d_alpha_dt, float d_beta_dt[3]) ;
		void combineStages(Cell2D &cell, float dt);
		void initialize_grid(int Nr, int Ntheta, float r_min, float r_max, float theta_min, float theta_max);
		float computeMaxSpeed();
		float computeCFL_dt(float CFL);
		void compute_constraints(Grid &grid_obj, int i, int j, int k, float &hamiltonian, float momentum[3]);
		void compute_time_derivatives(Grid &grid_obj, int i, int j, int k);
		void allocateGlobalGrid();
		void initializeData_Minkowski();
		void initializeKerrData(Grid &grid_obj);
		void initializeBinaryKerrData(Grid &grid_obj);
		float partialX_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);
		float partialY_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);
		float partialZ_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low);
		float computeTraceK(Grid &grid, int i, int j, int k);
		void compute_gauge_derivatives(Grid &grid_obj, int i, int j, int k, float &d_alpha_dt, float d_beta_dt[3]);
		void injectTTWave(Cell2D &cell, float x, float y, float z, float t);
		void solve_lichnerowicz(int max_iter, float tol, float dx, float dy, float dz);
		Cell2D& getCell(int i, int j, int k) {
			return globalGrid[i][j][k];
		}
		void export_Atildedt_slide(Grid &grid_obj, float time);
	private:
		std::vector<std::vector<std::vector<Grid::Cell2D>>> globalGrid;

};


float partialXX_alpha(Grid &grid_obj, int i, int j, int k);
float partialYY_alpha(Grid &grid_obj, int i, int j, int k);
float partialZZ_alpha(Grid &grid_obj, int i, int j, int k);
float second_partial_alpha(Grid &grid_obj, int i, int j, int k, int a, int b);
bool invert_3x3(const float m[3][3], float inv[3][3]);
void apply_boundary_conditions(Grid &grid_obj);
void export_K_3D(Grid &grid_obj);
void export_alpha_slice(Grid &grid_obj, int j);
void export_gauge_slice(Grid &grid_obj, int j);
void export_K_slice(Grid &grid_obj, int j);
void export_gamma_slice(Grid &grid_obj, int j, float time);
void export_tilde_gamma_3D(Grid &grid_obj);
void export_gauge_slice_anim(Grid &grid_obj, int j, float time);
void export_gauge_slice_xy(Grid &grid, int k, float time);
