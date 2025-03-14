#include <Geodesics.h>


void apply_boundary_conditions(Grid &grid_obj) {
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
            cell.gamma0[a][b] = cell.gamma[a][b];
            cell.K0[a][b]     = cell.K[a][b];
        }
    }
    cell.alpha0 = cell.alpha;
    for (int m = 0; m < 3; m++) {
        cell.beta0[m] = cell.beta[m];
    }
}

void Grid::updateIntermediateState(Cell2D &cell, double dtCoeff, int stageIndex) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.gamma[a][b] = cell.gamma0[a][b] + dtCoeff * cell.gammaStage[stageIndex][a][b];
            cell.K[a][b]     = cell.K0[a][b]     + dtCoeff * cell.KStage[stageIndex][a][b];
        }
    }
    cell.alpha = cell.alpha0 + dtCoeff * cell.alphaStage[stageIndex];
    for (int m = 0; m < 3; m++) {
        cell.beta[m] = cell.beta0[m] + dtCoeff * cell.betaStage[stageIndex][m];
    }
}

void Grid::storeStage(Cell2D &cell, int stage, double d_alpha_dt, double d_beta_dt[3]) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.gammaStage[stage][a][b] = cell.dgt[a][b];
            cell.KStage[stage][a][b]     = cell.dKt[a][b];
        }
    }
    cell.alphaStage[stage] = d_alpha_dt;
    for (int m = 0; m < 3; m++) {
        cell.betaStage[stage][m] = d_beta_dt[m];
    }
}

void Grid::combineStages(Cell2D &cell, double dt) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.gamma[a][b] = cell.gamma0[a][b] +
                (dt / 6.0) * (cell.gammaStage[0][a][b] +
                              2.0 * cell.gammaStage[1][a][b] +
                              2.0 * cell.gammaStage[2][a][b] +
                              cell.gammaStage[3][a][b]);
            cell.K[a][b] = cell.K0[a][b] +
                (dt / 6.0) * (cell.KStage[0][a][b] +
                              2.0 * cell.KStage[1][a][b] +
                              2.0 * cell.KStage[2][a][b] +
                              cell.KStage[3][a][b]);
        }
    }
    cell.alpha = cell.alpha0 +
        (dt / 6.0) * (cell.alphaStage[0] +
                      2.0 * cell.alphaStage[1] +
                      2.0 * cell.alphaStage[2] +
                      cell.alphaStage[3]);
    for (int m = 0; m < 3; m++) {
        cell.beta[m] = cell.beta0[m] +
            (dt / 6.0) * (cell.betaStage[0][m] +
                          2.0 * cell.betaStage[1][m] +
                          2.0 * cell.betaStage[2][m] +
                          cell.betaStage[3][m]);
    }
}
