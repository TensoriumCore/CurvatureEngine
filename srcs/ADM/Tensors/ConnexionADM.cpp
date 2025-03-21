#include "GridTensor.h"
#include <Geodesics.h>

/** 
 * This function compute Partial derivative of Tilde Gamma at each point of the grid
 * This is needed for BSSN implementations in the K_ij tensor and the 
 * Christoffel symbols
 * @param i , j , k the index of the cell
 * @param tildeGamma the output array of the tildeGamma tensor
 * @return void
 * */
void GridTensor::compute_dt_tildeGamma(Grid &grid_obj, int i, int j, int k, double dt_tildeGamma[3]) {
	Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
    
    double gammaInv[3][3], Atilde[3][3], dBeta[3][3];
    double alpha = cell.alpha;
    double beta[3] = { cell.beta[0], cell.beta[1], cell.beta[2] };

    int iP = std::min(i + 1, NX - 1);
    int iM = std::max(i - 1, 0);
    int jP = std::min(j + 1, NY - 1);
    int jM = std::max(j - 1, 0);
    int kP = std::min(k + 1, NZ - 1);
    int kM = std::max(k - 1, 0);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaInv[a][b] = cell.gamma_inv[a][b];
        }
    }

    double Ktrace = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * cell.K[a][b];
        }
    }

	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			cell.Atilde[a][b] = cell.K[a][b] - (1.0 / 3.0) * Ktrace * cell.tilde_gamma[a][b];
			cell.Atilde[a][b] *= cell.chi; 
		}
	}

	for (int n = 0; n < 3; n++) {
		dBeta[0][n] = (grid_obj.getCell(iP, j, k).beta[n] - grid_obj.getCell(iM, j, k).beta[n]) / (2.0 * DX);
		dBeta[1][n] = (grid_obj.getCell(i, jP, k).beta[n] - grid_obj.getCell(i, jM, k).beta[n]) / (2.0 * DY);
		dBeta[2][n] = (grid_obj.getCell(i, j, kP).beta[n] - grid_obj.getCell(i, j, kM).beta[n]) / (2.0 * DZ);
	}

 
	for (int i_comp = 0; i_comp < 3; i_comp++) {
		double div_Atilde = 0.0;
		double tildeGamma_Atilde = 0.0;
		double beta_term = 0.0;

		double KtraceP = 0.0, KtraceM = 0.0;
		for (int a = 0; a < 3; a++) {
			for (int b = 0; b < 3; b++) {
				KtraceP += grid_obj.getCell(iP, jP, kP).gamma_inv[a][b] * grid_obj.getCell(iP, jP, kP).K[a][b];
				KtraceM += grid_obj.getCell(iM, jM, kM).gamma_inv[a][b] * grid_obj.getCell(iM, jM, kM).K[a][b];
			}
		}

		for (int j_comp = 0; j_comp < 3; j_comp++) {
			double AtildeP = grid_obj.getCell(iP, j, k).K[i_comp][j_comp] -
				(1.0 / 3.0) * KtraceP * grid_obj.getCell(iP, j, k).tilde_gamma[i_comp][j_comp];

			double AtildeM = grid_obj.getCell(iM, j, k).K[i_comp][j_comp] -
				(1.0 / 3.0) * KtraceM * grid_obj.getCell(iM, j, k).tilde_gamma[i_comp][j_comp];

			div_Atilde += (AtildeP - AtildeM) / (2.0 * DX);
		}

		for (int j_comp = 0; j_comp < 3; j_comp++) {
			for (int k_comp = 0; k_comp < 3; k_comp++) {
				tildeGamma_Atilde += gammaInv[j_comp][k_comp] * Atilde[j_comp][k_comp];
			}
		}

		for (int j_comp = 0; j_comp < 3; j_comp++) {

			const auto &cellP = grid_obj.getCell(iP, j, k);
			const auto &cellM = grid_obj.getCell(iM, j, k);
			double GammaP = cellP.tildeGamma[i_comp];
			double GammaM = cellM.tildeGamma[i_comp];

			beta_term += grid_obj.getCell(i, j, k).beta[j_comp] * (GammaP - GammaM) / (2.0 * DX);
		}

		dt_tildeGamma[i_comp] = -2.0 * alpha * div_Atilde + 2.0 * alpha * tildeGamma_Atilde + beta_term;
	}
	for (int i_comp = 0; i_comp < 3; i_comp++) {
		dt_tildeGamma[i_comp] = cell.dt_tildeGamma[i_comp];
	}
}


/**
 * This function is a simmple implementation of the TildeGamma tensor
 * @param i , j , k the index of the cell
 * @param tildeGamma the output array of the tildeGamma tensor
 * @return void
 * */

void GridTensor::compute_tildeGamma(Grid &grid_obj, int i, int j, int k, double tildeGamma[3]) {
    double gammaInv[3][3];
    double christof[3][3][3];
    Grid::Cell2D &cell = grid_obj.getCell(i, j, k);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaInv[a][b] = cell.tildgamma_inv[a][b];
            for (int c = 0; c < 3; c++) {
                christof[a][b][c] = cell.Christoffel[a][b][c];
            }
        }
    }

    std::fill_n(tildeGamma, 3, 0.0);

	/*
	 * This use the christoffel symbols to compute the tildeGamma tensor
	 * \tilde{\Gamma}^i = \tilde_gamma^{jk} \Gamma^i_{jk}
	 * where \Gamma^i_{jk} is the christoffel symbol
	 * Then we contract the tensor conexion using the invert of the metric gamma_ij --> \tilde_gamma^{ij}
	 * */

    for (int i_comp = 0; i_comp < 3; i_comp++) { 
        for (int j_comp = 0; j_comp < 3; j_comp++) {
            for (int k_comp = 0; k_comp < 3; k_comp++) {
                tildeGamma[i_comp] += gammaInv[j_comp][k_comp] * christof[i_comp][j_comp][k_comp];
            }
        }
    }
	for (int i_comp = 0; i_comp < 3; i_comp++) {
		tildeGamma[i_comp] = cell.tildeGamma[i_comp];
	}
}

/** this function compute the Christoffel symbols at each point of the 3D grid
 * @param i , j , k the index of the cell
 * @param christof the output array of the christoffel tensor
 * @return void
 * */

void print_matrix_2D(const char *name, double matrix[3][3]) {
	printf("%s\n", name);
	for (int i = 0; i < 3; i++) {
		printf("[");
		for (int j = 0; j < 3; j++) {
			printf("%f ", matrix[i][j]);
		}
		printf("]\n");
	}
}

void GridTensor::compute_christoffel_3D(Grid &grid_obj, int i, int j, int k, double christof[3][3][3]) {
    Matrix matrix_obj;
    double g[3][3];
    double invg[3][3]; 

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            invg[a][b] = grid_obj.getCell(i, j, k).tildgamma_inv[a][b];
            g[a][b] = grid_obj.getCell(i, j, k).tilde_gamma[a][b];
            for (int c = 0; c < 3; c++) {
                christof[a][b][c] = grid_obj.getCell(i, j, k).Christoffel[a][b][c];
            }
        }
    }

    double dgamma[3][3][3];


	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			dgamma[0][a][b] = (i >= 2 && i <= NX - 3) ? fourth_order_diff(
					grid_obj.getCell(i+2, j, k).tilde_gamma[a][b],
					grid_obj.getCell(i+1, j, k).tilde_gamma[a][b],
					grid_obj.getCell(i-1, j, k).tilde_gamma[a][b],
					grid_obj.getCell(i-2, j, k).tilde_gamma[a][b],
					DX
					) : (i >= 1 && i <= NX - 2) ? second_order_diff(
						grid_obj.getCell(i+1, j, k).tilde_gamma[a][b],
						grid_obj.getCell(i-1, j, k).tilde_gamma[a][b],
						DX
						) : (i == 0) ? 
					(grid_obj.getCell(i+1, j, k).tilde_gamma[a][b] - grid_obj.getCell(i, j, k).tilde_gamma[a][b]) / DX
					: 
					(grid_obj.getCell(i, j, k).tilde_gamma[a][b] - grid_obj.getCell(i-1, j, k).tilde_gamma[a][b]) / DX;

			dgamma[1][a][b] = (j >= 2 && j <= NY - 3) ? fourth_order_diff(
					grid_obj.getCell(i, j+2, k).tilde_gamma[a][b],
					grid_obj.getCell(i, j+1, k).tilde_gamma[a][b],
					grid_obj.getCell(i, j-1, k).tilde_gamma[a][b],
					grid_obj.getCell(i, j-2, k).tilde_gamma[a][b],
					DY
					) : (j >= 1 && j <= NY - 2) ? second_order_diff(
						grid_obj.getCell(i, j+1, k).tilde_gamma[a][b],
						grid_obj.getCell(i, j-1, k).tilde_gamma[a][b],
						DY
						) : (j == 0) ? 
					(grid_obj.getCell(i, j+1, k).tilde_gamma[a][b] - grid_obj.getCell(i, j, k).tilde_gamma[a][b]) / DY
					: 
					(grid_obj.getCell(i, j, k).tilde_gamma[a][b] - grid_obj.getCell(i, j-1, k).tilde_gamma[a][b]) / DY;

			dgamma[2][a][b] = (k >= 2 && k <= NZ - 3) ? fourth_order_diff(
					grid_obj.getCell(i, j, k+2).tilde_gamma[a][b],
					grid_obj.getCell(i, j, k+1).tilde_gamma[a][b],
					grid_obj.getCell(i, j, k-1).tilde_gamma[a][b],
					grid_obj.getCell(i, j, k-2).tilde_gamma[a][b],
					DZ
					) : (k >= 1 && k <= NZ - 2) ? second_order_diff(
						grid_obj.getCell(i, j, k+1).tilde_gamma[a][b],
						grid_obj.getCell(i, j, k-1).tilde_gamma[a][b],
						DZ
						) : (k == 0) ? 
					(grid_obj.getCell(i, j, k+1).tilde_gamma[a][b] - grid_obj.getCell(i, j, k).tilde_gamma[a][b]) / DZ
					: 
					(grid_obj.getCell(i, j, k).tilde_gamma[a][b] - grid_obj.getCell(i, j, k-1).tilde_gamma[a][b]) / DZ;

			/* if (i == NX/2 && j == NY/2 && k == NZ/2) { */
			/* 	printf("∆X = %e, ∆Y = %e, ∆Z = %e\n", DX, DY, DZ); */
			/* 	printf("∂_x g[0][0] = %e\n", dgamma[0][0][0]); */
			/* 	printf("∂_y g[0][0] = %e\n", dgamma[1][0][0]); */
			/* 	printf("∂_z g[0][0] = %e\n", dgamma[2][0][0]); */
			/* } */


		}
	}

    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
                double sum = 0.0;
                for (int ll = 0; ll < 3; ll++) {
                    double tmp = dgamma[aa][ll][bb] + dgamma[bb][ll][aa] - dgamma[ll][aa][bb];
                    sum += invg[kk][ll] * tmp;
                }
                grid_obj.getCell(i, j, k).Christoffel[kk][aa][bb] = 0.5 * sum;
            }
        }
    }
}
