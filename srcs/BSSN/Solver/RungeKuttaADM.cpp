#include <Geodesics.h>

void Grid::logger_evolve(Grid &grid_obj, double dt, int nstep)
{    
	Cell2D &cell = grid_obj.getCell(0, 0, 0);
	printf("=============================================\n");
	printf("Time step: %f, nstep: %d\n", dt, nstep);
	printf("=============================================\n");
	printf("alpha: %e\n", cell.gauge.alpha);
	printf("beta: %e %e %e\n", cell.gauge.beta[0], cell.gauge.beta[1], cell.gauge.beta[2]);
	printf("Hamiltonian: %e\n", cell.matter.hamiltonian);
	printf("Momentum: %e %e %e\n", cell.matter.momentum[0], cell.matter.momentum[1], cell.matter.momentum[2]);
	printf("dtAtilde:\n");
	for (int i = 0; i < 3; i++)
		printf("	  %e %e %e\n", cell.atilde.dt_Atilde[i][0], \
								   cell.atilde.dt_Atilde[i][1], 
								   cell.atilde.dt_Atilde[i][2]);
	printf("Chi: %e\n", cell.dt_chi);
	printf("=============================================\n");
	printf("\n\n");
}



void Grid::evolve(Grid &grid_obj, double dtInitial, int nSteps) {
    initialize_grid();
    GridTensor gridTensor;
    double CFL = 0.5;
    double dt = dtInitial;
    double hamiltonian;
    double momentum[3];
	double L = 9.0;
    double x_min = -L, x_max = L;
    double y_min = -L, y_max = L;
    double z_min = -L, z_max = L;


    double dx = (x_max - x_min) / (NX - 1);
    double dy = (y_max - y_min) / (NY - 1);
    double dz = (z_max - z_min) / (NZ - 1);
    for (int step = 0; step < nSteps; step++) {
        dt = computeCFL_dt(CFL);
        apply_boundary_conditions(grid_obj);

#pragma omp parallel
        {
            auto forEachCell = [&](auto func) {
#pragma omp for collapse(3) schedule(runtime)
                for (int i = 1; i < NX - 1; i++) {
                    for (int j = 1; j < NY - 1; j++) {
                        for (int k = 1; k < NZ - 1; k++) {
							double x = x_min + i*DX;
                            double y = y_min + j*DY;
                            double z = z_min + k*DZ;
                            injectTTWave(globalGrid[i][j][k], x, y, z, grid_obj.time);
                            func(i, j, k);
                        }
                    }
                }
            };

            forEachCell([&](int i, int j, int k) {
                copyInitialState(globalGrid[i][j][k]);
            });

            forEachCell([&](int i, int j, int k) {
                compute_time_derivatives(grid_obj, i, j, k);
                double d_alpha_dt, d_beta_dt[3];
                compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);
                compute_constraints(grid_obj, i, j, k, hamiltonian, momentum);
                storeStage(globalGrid[i][j][k], 0, d_alpha_dt, d_beta_dt);
            });

            forEachCell([&](int i, int j, int k) {
                updateIntermediateState(globalGrid[i][j][k], 0.5 * dt, 0);
            });

            forEachCell([&](int i, int j, int k) {
                compute_time_derivatives(grid_obj, i, j, k);
                double d_alpha_dt, d_beta_dt[3];
                compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);
                storeStage(globalGrid[i][j][k], 1, d_alpha_dt, d_beta_dt);
            });

            forEachCell([&](int i, int j, int k) {
                updateIntermediateState(globalGrid[i][j][k], 0.5 * dt, 1);
            });

            forEachCell([&](int i, int j, int k) {
                compute_time_derivatives(grid_obj, i, j, k);
                double d_alpha_dt, d_beta_dt[3];
                compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);
                compute_constraints(grid_obj, i, j, k, hamiltonian, momentum);
                storeStage(globalGrid[i][j][k], 2, d_alpha_dt, d_beta_dt);
            });

            forEachCell([&](int i, int j, int k) {
                updateIntermediateState(globalGrid[i][j][k], dt, 2);
            });

            forEachCell([&](int i, int j, int k) {
                compute_time_derivatives(grid_obj, i, j, k);
                double d_alpha_dt, d_beta_dt[3];
                compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);
                compute_constraints(grid_obj, i, j, k, hamiltonian, momentum);
                storeStage(globalGrid[i][j][k], 3, d_alpha_dt, d_beta_dt);
            });

            forEachCell([&](int i, int j, int k) {
                combineStages(globalGrid[i][j][k], dt);
            });
        }

#pragma omp single nowait
        {
            logger_evolve(grid_obj, dt, step);
			double current_time = step * dt;
			export_gamma_slice(grid_obj, NY / 2, dt);
			grid_obj.appendConstraintL2ToCSV("constraints_evolution.csv", current_time);
			if (step == nSteps - 1) {
                printf("Exporting slices\n");
                export_K_slice(grid_obj, NY / 2);
                export_gauge_slice(grid_obj, NY / 2);
                gridTensor.export_christoffel_slice(grid_obj, NX / 2);
                export_K_3D(grid_obj);
            }
        }

        grid_obj.time += dt;
    }
}
