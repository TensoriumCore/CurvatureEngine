#include <Geodesics.h>

/** 
 * This function compute Partial derivative of Tilde Gamma at each point of the grid
 * This is needed for BSSN implementations in the K_ij tensor and the 
 * conn.Christoffel symbols
 * @param i , j , k the index of the cell
 * @param tildeGamma the output array of the tildeGamma tensor
 * @return void
 * */
void GridTensor::compute_dt_tildeGamma(Grid &grid_obj, int i, int j, int k, double dt_tildeGamma[3]) {
	Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
    
    double gammaInv[3][3], Atilde[3][3], dBeta[3][3];
    double alpha = cell.gauge.alpha;
    double beta[3] = { cell.gauge.beta[0], cell.gauge.beta[1], cell.gauge.beta[2] };

    int iP = std::min(i + 1, NX - 1);
    int iM = std::max(i - 1, 0);
    int jP = std::min(j + 1, NY - 1);
    int jM = std::max(j - 1, 0);
    int kP = std::min(k + 1, NZ - 1);
    int kM = std::max(k - 1, 0);


	gammaInv[0][0] = cell.geom.tildgamma_inv[0][0];
	gammaInv[0][1] = cell.geom.tildgamma_inv[0][1];
	gammaInv[0][2] = cell.geom.tildgamma_inv[0][2];
	gammaInv[1][0] = cell.geom.tildgamma_inv[1][0];
	gammaInv[1][1] = cell.geom.tildgamma_inv[1][1];
	gammaInv[1][2] = cell.geom.tildgamma_inv[1][2];
	gammaInv[2][0] = cell.geom.tildgamma_inv[2][0];
	gammaInv[2][1] = cell.geom.tildgamma_inv[2][1];
	gammaInv[2][2] = cell.geom.tildgamma_inv[2][2];

    double Ktrace = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * cell.curv.K[a][b];
        }
    }

	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			cell.atilde.Atilde[a][b] = cell.curv.K[a][b] - (1.0 / 3.0) * Ktrace * cell.geom.tilde_gamma[a][b];
			cell.atilde.Atilde[a][b] *= cell.chi; 
		}
	}

	for (int n = 0; n < 3; n++) {
		dBeta[0][n] = (grid_obj.getCell(iP, j, k).gauge.beta[n] - grid_obj.getCell(iM, j, k).gauge.beta[n]) / (2.0 * DX);
		dBeta[1][n] = (grid_obj.getCell(i, jP, k).gauge.beta[n] - grid_obj.getCell(i, jM, k).gauge.beta[n]) / (2.0 * DY);
		dBeta[2][n] = (grid_obj.getCell(i, j, kP).gauge.beta[n] - grid_obj.getCell(i, j, kM).gauge.beta[n]) / (2.0 * DZ);
	}

 
	for (int i_comp = 0; i_comp < 3; i_comp++) {
		double div_Atilde = 0.0;
		double tildeGamma_Atilde = 0.0;
		double beta_term = 0.0;

		double KtraceP = 0.0, KtraceM = 0.0;
		for (int a = 0; a < 3; a++) {
			for (int b = 0; b < 3; b++) {
				KtraceP += grid_obj.getCell(iP, jP, kP).geom.gamma_inv[a][b] * grid_obj.getCell(iP, jP, kP).curv.K[a][b];
				KtraceM += grid_obj.getCell(iM, jM, kM).geom.gamma_inv[a][b] * grid_obj.getCell(iM, jM, kM).curv.K[a][b];
			}
		}

		for (int j_comp = 0; j_comp < 3; j_comp++) {
			double AtildeP = grid_obj.getCell(iP, j, k).curv.K[i_comp][j_comp] -
				(1.0 / 3.0) * KtraceP * grid_obj.getCell(iP, j, k).geom.tilde_gamma[i_comp][j_comp];

			double AtildeM = grid_obj.getCell(iM, j, k).curv.K[i_comp][j_comp] -
				(1.0 / 3.0) * KtraceM * grid_obj.getCell(iM, j, k).geom.tilde_gamma[i_comp][j_comp];

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
			double GammaP = cell.conn.tildeGamma[i_comp];
			double GammaM = cell.conn.tildeGamma[i_comp];

			beta_term += grid_obj.getCell(i, j, k).gauge.beta[j_comp] * (GammaP - GammaM) / (2.0 * DX);
		}

		dt_tildeGamma[i_comp] = -2.0 * alpha * div_Atilde + 2.0 * alpha * tildeGamma_Atilde + beta_term;
	}
	for (int i_comp = 0; i_comp < 3; i_comp++) {
		dt_tildeGamma[i_comp] = cell.geom.dt_tildeGamma[i_comp];
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

    /* for (int a = 0; a < 3; a++) { */
    /*     for (int b = 0; b < 3; b++) { */
    /*         gammaInv[a][b] = cell.geom.tildgamma_inv[a][b]; */
    /*         for (int c = 0; c < 3; c++) { */
    /*             christof[a][b][c] = cell.conn.Christoffel[a][b][c]; */
    /*         } */
    /*     } */
    /* } */

    std::fill_n(tildeGamma, 3, 0.0);

	/*
	 * This use the christoffel symbols to compute the tildeGamma tensor
	 * \tilde{\Gamma}^i = \tilde_gamma^{jk} \Gamma^i_{jk}
	 * where \Gamma^i_{jk} is the christoffel symbol
	 * Then we contract the tensor conexion using the invert of the metric gamma_ij --> \tilde_gamma^{ij}
	 * */
#pragma omp simd collapse(3)
    for (int i_comp = 0; i_comp < 3; i_comp++) { 
        for (int j_comp = 0; j_comp < 3; j_comp++) {
            for (int k_comp = 0; k_comp < 3; k_comp++) {
                tildeGamma[i_comp] += cell.geom.tildgamma_inv[j_comp][k_comp] * cell.conn.Christoffel[i_comp][j_comp][k_comp];
            }
        }
    }
    cell.conn.tildeGamma[0] = tildeGamma[0];
    cell.conn.tildeGamma[1] = tildeGamma[1];
    cell.conn.tildeGamma[2] = tildeGamma[2];
}

/** this function compute the conn.Christoffel symbols at each point of the 3D grid
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
    const auto& cell = grid_obj.getCell(i, j, k);

    double dgamma[3][3][3];
    auto get_g = [&](int ii, int jj, int kk, int a, int b) -> double {
        return grid_obj.getCell(ii, jj, kk).geom.tilde_gamma[a][b];
    };

    bool in_x = (i >= 2 && i <= NX - 3);
    bool in_y = (j >= 2 && j <= NY - 3);
    bool in_z = (k >= 2 && k <= NZ - 3);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            if (in_x) {
                dgamma[0][a][b] = fourth_order_diff(
                    get_g(i+2,j,k,a,b), get_g(i+1,j,k,a,b),
                    get_g(i-1,j,k,a,b), get_g(i-2,j,k,a,b),
                    DX
                );
            } else if (i >= 1 && i <= NX - 2) {
                dgamma[0][a][b] = second_order_diff(
                    get_g(i+1,j,k,a,b), get_g(i-1,j,k,a,b),
                    DX
                );
            } else if (i == 0) {
				dgamma[0][a][b] = second_order_diff(
					get_g(i+1,j,k,a,b), get_g(i,j,k,a,b),
					DX
				);
            } else {
                dgamma[0][a][b] = second_order_diff(
					get_g(i,j,k,a,b), get_g(i-1,j,k,a,b),
					DX
				);
            }

            // ∂_y
            if (in_y) {
                dgamma[1][a][b] = fourth_order_diff(
                    get_g(i,j+2,k,a,b), get_g(i,j+1,k,a,b),
                    get_g(i,j-1,k,a,b), get_g(i,j-2,k,a,b),
                    DY
                );
            } else if (j >= 1 && j <= NY - 2) {
                dgamma[1][a][b] = second_order_diff(
                    get_g(i,j+1,k,a,b), get_g(i,j-1,k,a,b),
                    DY
                );
            } else if (j == 0) {
                dgamma[1][a][b] = second_order_diff(
					get_g(i,j+1,k,a,b), get_g(i,j,k,a,b),
					DY
				); 
            } else {
                dgamma[1][a][b] = second_order_diff(
					get_g(i,j,k,a,b), get_g(i,j-1,k,a,b),
					DY
				);
			}
            if (in_z) {
                dgamma[2][a][b] = fourth_order_diff(
                    get_g(i,j,k+2,a,b), get_g(i,j,k+1,a,b),
                    get_g(i,j,k-1,a,b), get_g(i,j,k-2,a,b),
                    DZ
                );
            } else if (k >= 1 && k <= NZ - 2) {
                dgamma[2][a][b] = second_order_diff(
                    get_g(i,j,k+1,a,b), get_g(i,j,k-1,a,b),
                    DZ
                );
            } else if (k == 0) {
                dgamma[2][a][b] = second_order_diff( 
					get_g(i,j,k+1,a,b), get_g(i,j,k,a,b),
					DZ
				);
            } else {
                dgamma[2][a][b] = second_order_diff(
					get_g(i,j,k,a,b), get_g(i,j,k-1,a,b), 
					DZ
				);
			}
        }
    }
#pragma omp simd collapse(3)
    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
                double sum = 0.0;
                for (int ll = 0; ll < 3; ll++) {
                    double tmp = dgamma[aa][ll][bb] + dgamma[bb][ll][aa] - dgamma[ll][aa][bb];
                    sum += cell.geom.tildgamma_inv[kk][ll] * tmp;
                }
                christof[kk][aa][bb] = 0.5 * sum;
                grid_obj.getCell(i, j, k).conn.Christoffel[kk][aa][bb] = christof[kk][aa][bb];
            }
        }
    }
}
