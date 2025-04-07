#include <Geodesics.h>

void Grid::compute_gauge_derivatives(Grid &grid_obj, int i, int j, int k, float &d_alpha_dt, float d_beta_dt[3]) {
    Grid::Cell2D &cell = globalGrid[i][j][k];
    float gammaLocal[3][3], KLocal[3][3];
	GridTensor gridTensor;
    float gammaInv[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaLocal[a][b] = cell.geom.gamma[a][b];
            KLocal[a][b]     = cell.curv.K[a][b];
			gammaInv[a][b] = cell.geom.gamma_inv[a][b];	
        }
    }

    float Ktrace = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * KLocal[a][b];
		}
    }

    d_alpha_dt = -2.0 * cell.gauge.alpha * Ktrace ;

    float eta = 2.0 / (1.0 + std::fabs(Ktrace));
    float d_Gamma_dt[3] = {0.0, 0.0, 0.0}; 

    float tildeGamma[3];
    gridTensor.compute_tildeGamma(grid_obj, i, j, k, tildeGamma); 
    gridTensor.compute_dt_tildeGamma(grid_obj, i, j, k, d_Gamma_dt);

    for (int m = 0; m < 3; m++) {
        d_beta_dt[m] = 3.0 / 4.0 * d_Gamma_dt[m] - eta * cell.gauge.beta[m];
	}
	cell.gauge.dt_alpha = d_alpha_dt;
	for (int m = 0; m < 3; m++) {
		cell.gauge.dt_beta[m] = d_beta_dt[m];
	}
}

