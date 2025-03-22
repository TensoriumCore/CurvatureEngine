#include <Geodesics.h>

void Grid::logger_evolve(Grid &grid_obj, double dt, int nstep)
{    
	Cell2D &cell = grid_obj.getCell(0, 0, 0);
	printf("=====================================\n");
	printf("Time step: %f, nstep: %d\n", dt, nstep);
	printf("=====================================\n");
	printf("alpha: %f\n", cell.alpha);
	printf("beta: %f %f %f\n", cell.beta[0], cell.beta[1], cell.beta[2]);
	printf("Hamiltonian: %f\n", cell.hamiltonian);
	printf("Momentum: %f %f %f\n", cell.momentum[0], cell.momentum[1], cell.momentum[2]);
	printf("Ricciscalar: %f\n", cell.Ricci);
	printf("=====================================\n");
	printf("\n\n");
	
}

void Grid::evolve(Grid &grid_obj, double dtinitital, int nSteps) {
	initialize_grid(); 
	GridTensor gridTensor;
	double CFL = 0.5;
	double dt = dtinitital;
	double hamiltonian, momentum[3];

	for (int step = 0; step < nSteps; step++) {
		dt = computeCFL_dt(CFL);
		apply_boundary_conditions(grid_obj);
#pragma omp parallel 
		{
#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						copyInitialState(globalGrid[i][j][k]);
					}
				}
			}
#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						compute_time_derivatives(grid_obj, i, j, k);
						double d_alpha_dt, d_beta_dt[3];
						compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);
						compute_constraints(grid_obj, i, j, k, hamiltonian, momentum);
						storeStage(globalGrid[i][j][k], 0, d_alpha_dt, d_beta_dt);
					}
				}
			}

#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						updateIntermediateState(globalGrid[i][j][k], 0.5 * dt, 0);
					}
				}
			}
#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						compute_time_derivatives(grid_obj, i, j, k);
						double d_alpha_dt, d_beta_dt[3];
						compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);
						storeStage(globalGrid[i][j][k], 1, d_alpha_dt, d_beta_dt);
					}
				}
			}
#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						updateIntermediateState(globalGrid[i][j][k], 0.5 * dt, 1);
					}
				}
			}
#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						compute_time_derivatives(grid_obj, i, j, k);
						double d_alpha_dt, d_beta_dt[3];
						compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);
						compute_constraints(grid_obj, i, j, k, hamiltonian, momentum);
						storeStage(globalGrid[i][j][k], 2, d_alpha_dt, d_beta_dt);
					}
				}
			}
#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						updateIntermediateState(globalGrid[i][j][k], dt, 2);
					}
				}
			}
#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						compute_time_derivatives(grid_obj, i, j, k);
						double d_alpha_dt, d_beta_dt[3];
						compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);
						compute_constraints(grid_obj, i, j, k, hamiltonian, momentum);
						storeStage(globalGrid[i][j][k], 3, d_alpha_dt, d_beta_dt);
					}
				}
			}
#pragma omp for collapse(3) schedule(dynamic)
			for (int i = 1; i < NX - 1; i++) {
				for (int j = 1; j < NY - 1; j++) {
					for (int k = 1; k < NZ - 1; k++) {
						combineStages(globalGrid[i][j][k], dt);
					}
				}
			}
#pragma omp single
			{
				logger_evolve(grid_obj, dt, step);
				export_chi_slice(grid_obj, step * dt);
				export_gamma_slice(grid_obj, NY / 2, step * dt);
				if (step == nSteps- 1)
				{
					printf("Exporting slices\n");
					export_K_slice(grid_obj, NY / 2);
					export_gauge_slice(grid_obj, NY / 2);
					gridTensor.export_christoffel_slice(grid_obj, NY / 2);
					export_alpha_slice(grid_obj, NY / 2);
					export_K_3D(grid_obj);
					export_tilde_gamma_3D(grid_obj);
					/* export_constraints("Output/constraints_N128.csv"); */
				}
			}
		}
	}
}
