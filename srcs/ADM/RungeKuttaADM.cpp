#include <Geodesics.h>



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
						compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);
						/* update_fluid_velocity(i, j, k, dt); */
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
						compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);
						/* update_fluid_velocity(i, j, k, dt); */
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
						compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);
						/* update_fluid_velocity(i, j, k, dt); */
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
						compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);
						/* update_fluid_velocity(i, j, k, dt); */
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
				if (step == nSteps- 1)
				{
					printf("Exporting slices\n");
					export_K_slice(grid_obj, NY / 2);
					export_gauge_slice(grid_obj, NY / 2);
					export_gamma_slice(grid_obj, NY / 2);
					gridTensor.export_christoffel_slice(grid_obj, NY / 2);
					/* export_fluid_slice(NY / 2); */
					/* export_energy_momentum_tensor_slice(NY / 2); */
					export_K_3D(grid_obj);
				}
			}
		}
	}
}
