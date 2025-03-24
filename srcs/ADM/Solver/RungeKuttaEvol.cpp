#include <Geodesics.h>


void apply_asymptotic_boundary_conditions(Grid &grid_obj, double R_max) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = i * DX ;
                double y = j * DY ;
                double z = k * DZ ;
                double r = sqrt(x * x + y * y + z * z);

                if (r > R_max) {
                    grid_obj.getCell(i, j, k).gauge.alpha = 1.0;
                    grid_obj.getCell(i, j, k).gauge.beta[0] = 0.0;
                    grid_obj.getCell(i, j, k).gauge.beta[1] = 0.0;
                    grid_obj.getCell(i, j, k).gauge.beta[2] = 0.0;

                    for (int a = 0; a < 3; a++) {
                        for (int b = 0; b < 3; b++) {
                            grid_obj.getCell(i, j, k).curv.K[a][b] = 0.0;
                        }
                    }
                }
            }
        }
    }
}

/* void apply_horizon_excision(Grid &grid_obj, double r_H) { */
/*     for (int i = 0; i < NX; i++) { */
/*         for (int j = 0; j < NY; j++) { */
/*             for (int k = 0; k < NZ; k++) { */
/*                 double x = i * DX; */
/*                 double y = j * DY; */
/*                 double z = k * DZ; */
/*                 double r = sqrt(x * x + y * y + z * z); */
/*  */
/*                 if (r < r_H) { */
/*                     grid_obj.getCell(i, j, k).gauge.alpha = 1.0 / sqrt(1.0 + 2.0 * grid_obj.getCell(i, j, k).H); */
/*                     grid_obj.getCell(i, j, k).gauge.beta[0] = grid_obj.getCell(i, j, k).l_x; */
/*                     grid_obj.getCell(i, j, k).gauge.beta[1] = grid_obj.getCell(i, j, k).l_y; */
/*                     grid_obj.getCell(i, j, k).gauge.beta[2] = grid_obj.getCell(i, j, k).l_z; */
/*                 } */
/*             } */
/*         } */
/*     } */
/* } */
/*  */

void apply_boundary_conditions(Grid &grid_obj) {
	apply_asymptotic_boundary_conditions(grid_obj, 128.0);
    for (int j = 0; j < NY; j++) {
        for (int k = 0; k < NZ; k++) {
            grid_obj.getCell(0, j, k) = grid_obj.getCell(1, j, k);
            grid_obj.getCell(NX - 1, j, k) = grid_obj.getCell(NX - 2, j, k);
        }
    }

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            grid_obj.getCell(i, 0, k) = grid_obj.getCell(i, 1, k);
            grid_obj.getCell(i, NY - 1, k) = grid_obj.getCell(i, NY - 2, k);
        }
    }

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            grid_obj.getCell(i, j, 0) = grid_obj.getCell(i, j, 1);
            grid_obj.getCell(i, j, NZ - 1) = grid_obj.getCell(i, j, NZ - 2);
        }
    }
}


void Grid::copyInitialState(Cell2D &cell) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.geom.tilde_gamma0[a][b] = cell.geom.gamma[a][b];
            cell.curv.K0[a][b]     = cell.curv.K[a][b];
			cell.atilde.Atilde0[a][b] = cell.atilde.Atilde[a][b];
        }
    }
    cell.gauge.alpha0 = cell.gauge.alpha;
    for (int m = 0; m < 3; m++) {
        cell.gauge.beta0[m] = cell.gauge.beta[m];
    }
}

void Grid::updateIntermediateState(Cell2D &cell, double dtCoeff, int stageIndex) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.geom.tilde_gamma[a][b] = cell.geom.tilde_gamma0[a][b] + dtCoeff * cell.gammaStage[stageIndex][a][b];
            cell.curv.K[a][b]     = cell.curv.K0[a][b]     + dtCoeff * cell.KStage[stageIndex][a][b];
        }
    }
    cell.gauge.alpha = cell.gauge.alpha0 + dtCoeff * cell.gauge.alphaStage[stageIndex];
    for (int m = 0; m < 3; m++) {
        cell.gauge.beta[m] = cell.gauge.beta0[m] + dtCoeff * cell.gauge.betaStage[stageIndex][m];
    }
}

void Grid::storeStage(Cell2D &cell, int stage, double d_alpha_dt, double d_beta_dt[3]) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.gammaStage[stage][a][b] = cell.dgt[a][b];
            cell.KStage[stage][a][b]     = cell.curv.dKt[a][b];
        }
    }
    cell.gauge.alphaStage[stage] = d_alpha_dt;
    for (int m = 0; m < 3; m++) {
        cell.gauge.betaStage[stage][m] = d_beta_dt[m];
    }
}

void Grid::combineStages(Cell2D &cell, double dt) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.geom.tilde_gamma[a][b] = cell.geom.tilde_gamma0[a][b] +
                (dt / 6.0) * (cell.gammaStage[0][a][b] +
                              2.0 * cell.gammaStage[1][a][b] +
                              2.0 * cell.gammaStage[2][a][b] +
                              cell.gammaStage[3][a][b]);
            cell.curv.K[a][b] = cell.curv.K0[a][b] +
                (dt / 6.0) * (cell.KStage[0][a][b] +
                              2.0 * cell.KStage[1][a][b] +
                              2.0 * cell.KStage[2][a][b] +
                              cell.KStage[3][a][b]);

			cell.atilde.Atilde[a][b] = cell.atilde.Atilde0[a][b] +
				(dt / 6.0) * (cell.atilde.AtildeStage[0][a][b] +
						2.0 * cell.atilde.AtildeStage[1][a][b] +
						2.0 * cell.atilde.AtildeStage[2][a][b] +
						cell.atilde.AtildeStage[3][a][b]);

        }
    }
    cell.gauge.alpha = cell.gauge.alpha0 +
        (dt / 6.0) * (cell.gauge.alphaStage[0] +
                      2.0 * cell.gauge.alphaStage[1] +
                      2.0 * cell.gauge.alphaStage[2] +
                      cell.gauge.alphaStage[3]);
    for (int m = 0; m < 3; m++) {
        cell.gauge.beta[m] = cell.gauge.beta0[m] +
            (dt / 6.0) * (cell.gauge.betaStage[0][m] +
                          2.0 * cell.gauge.betaStage[1][m] +
                          2.0 * cell.gauge.betaStage[2][m] +
                          cell.gauge.betaStage[3][m]);
    }
}
