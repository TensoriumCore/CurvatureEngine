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
void GridTensor::compute_dt_tildeGamma(int i, int j, int k, double dt_tildeGamma[3]) {
	Grid grid_obj;
	Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
    
    double gammaInv[3][3], KLocal[3][3], Atilde[3][3], dBeta[3][3];
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
            KLocal[a][b] = cell.K[a][b];
        }
    }

    double Ktrace = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * KLocal[a][b];
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Atilde[a][b] = KLocal[a][b] - (1.0 / 3.0) * Ktrace * cell.gamma[a][b];
        }
    }

    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            dBeta[m][n] = (grid_obj.getCell(iP, jP, kP).beta[n] - grid_obj.getCell(iM, jM, kM).beta[n]) / (2.0 * DX);
        }
    }

    for (int i_comp = 0; i_comp < 3; i_comp++) {
        double div_Atilde = 0.0;
        double tildeGamma_Atilde = 0.0;
        double beta_term = 0.0;
		/*
		 * The divergence of the Atilde tensor is computed using the formula:
		 * \partial_t Atilde^i = \partial_j K^{ij}
		 * where K^{ij} is the extrinsic curvature tensor
		 */
        for (int j_comp = 0; j_comp < 3; j_comp++) {
            div_Atilde += (grid_obj.getCell(iP, jP, kP).K[i_comp][j_comp] -
                           grid_obj.getCell(iM, jM, kM).K[i_comp][j_comp]) / (2.0 * DX);
        }

        for (int j_comp = 0; j_comp < 3; j_comp++) {
            for (int k_comp = 0; k_comp < 3; k_comp++) {
                tildeGamma_Atilde += gammaInv[j_comp][k_comp] * Atilde[j_comp][k_comp];
            }
        }

        for (int j_comp = 0; j_comp < 3; j_comp++) {
            beta_term += beta[j_comp] * (grid_obj.getCell(iP, jP, kP).beta[i_comp] -
                                          grid_obj.getCell(iM, jM, kM).beta[i_comp]) / (2.0 * DX);
        }
		/*
		 * The tildeGamma tensor is computed using the formula:
		 * \partial_t \tilde{\Gamma}^i = -2 \alpha \partial_j K^{ij} + 2 \alpha \tilde{\Gamma}^i_{jk} K^{jk} + \beta^j \partial_j \tilde{\Gamma}^i
		 * using the gauge source function \beta^j = \partial_t \beta^j and the lapse function \alpha = \partial_t \alpha
		 */
        dt_tildeGamma[i_comp] = -2.0 * alpha * div_Atilde + 2.0 * alpha * tildeGamma_Atilde + beta_term;
    }
}


/**
 * This function is a simmple implementation of the TildeGamma tensor
 * @param i , j , k the index of the cell
 * @param tildeGamma the output array of the tildeGamma tensor
 * @return void
 * */

void GridTensor::compute_tildeGamma(int i, int j, int k, double tildeGamma[3]) {
    double gammaInv[3][3];
    double christof[3][3][3];
	Grid grid_obj;
    Grid::Cell2D &cell = grid_obj.getCell(i, j, k);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaInv[a][b] = cell.gamma_inv[a][b];
            for (int c = 0; c < 3; c++) {
                christof[a][b][c] = cell.Christoffel[a][b][c];
            }
        }
    }

    std::fill_n(tildeGamma, 3, 0.0);

	/*
	 * This use the christoffel symbols to compute the tildeGamma tensor
	 * \tilde{\Gamma}^i = \gamma^{jk} \Gamma^i_{jk}
	 * where \Gamma^i_{jk} is the christoffel symbol
	 * Then we contract the tensor conexion using the invert of the metric gamma_ij --> \gamma^{ij}
	 * */

    for (int i_comp = 0; i_comp < 3; i_comp++) { 
        for (int j_comp = 0; j_comp < 3; j_comp++) {
            for (int k_comp = 0; k_comp < 3; k_comp++) {
                tildeGamma[i_comp] += gammaInv[j_comp][k_comp] * christof[i_comp][j_comp][k_comp];
            }
        }
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
			invg[a][b] = grid_obj.getCell(i, j, k).gamma_inv[a][b];
			g[a][b] = grid_obj.getCell(i, j, k).gamma[a][b];
			for (int c = 0; c < 3; c++) {
				christof[a][b][c] = grid_obj.getCell(i, j, k).Christoffel[a][b][c];
			}
		}
	}
	
	/* if (i == NX/2 && j == NY/2 && k == NZ/2) { */
	/* 	print_matrix_2D("gamma", g); */
	/* 	print_matrix_2D("gamma_inv", invg); */
	/* } */
	/* printf("Vérification de g_ij * g^jk :\n"); */
	/* for (int i = 0; i < 3; i++) { */
	/* 	for (int k = 0; k < 3; k++) { */
	/* 		double sum = 0.0; */
	/* 		for (int j = 0; j < 3; j++) { */
	/* 			sum += grid_obj.getCell(i, j, k).gamma[i][j] * grid_obj.getCell(i, j, k).gamma_inv[j][k]; */
	/* 		} */
	/* 		printf("%f ", sum); */
	/* 	} */
	/* 	printf("\n"); */
	/* } */

	/*
	 * First we compute the partial derivative of the metric tensor
	 * Then we use the partial derivative of the metric tensor to compute the Christoffel symbols
	 * The Christoffel symbols are computed using the formula:
	 * \Gamma^i_{jk} = 1/2 g^{il} ( \partial_j g_{lk} + \partial_k g_{jl} - \partial_l g_{jk} )
	 * where g_{ij} is the metric tensor and g^{ij} is the inverse of the metric tensor
	 * */

	double dgamma[3][3][3];
	for(int a=0;a<3;a++){
		for(int b=0;b<3;b++){
			dgamma[0][a][b] = partialX_gamma(grid_obj, i, j, k, a, b);
			dgamma[1][a][b] = partialY_gamma(grid_obj, i, j, k, a, b);
			dgamma[2][a][b] = partialZ_gamma(grid_obj, i, j, k, a, b);
			/* printf("dx_g[%d][%d] = %e, dy_g[%d][%d] = %e, dz_g[%d][%d] = %e\n",  */
			/* 		a, b, partialX_gamma(grid_obj, i, j, k, a, b), */
			/* 		a, b, partialY_gamma(grid_obj, i, j, k, a, b), */
			/* 		a, b, partialZ_gamma(grid_obj, i, j, k, a, b)); */
			/* printf("Test cohérence : ∂xg[%d][%d] - ∂yg[%d][%d] = %e\n",  */
			/* 		a, b, a, b,  */
			/* 		partialX_gamma(grid_obj, i, j, k, a, b) - partialY_gamma(grid_obj, i, j, k, a, b)); */
			/*  */
			/* double avg_dx = (partialX_gamma(grid_obj, i, j, k, a, b) + partialX_gamma(grid_obj, i+1, j, k, a, b) + partialX_gamma(grid_obj, i-1, j, k, a, b)) / 3.0; */
			/* double avg_dy = (partialY_gamma(grid_obj, i, j, k, a, b) + partialY_gamma(grid_obj, i, j+1, k, a, b) + partialY_gamma(grid_obj, i, j-1, k, a, b)) / 3.0; */
			/* printf("Moyenne dx_g[%d][%d] = %e, Moyenne dy_g[%d][%d] = %e\n", a, b, avg_dx, a, b, avg_dy); */
		}
	}


	/*
	 * The Christoffel symbols are computed using the formula:
	 * \Gamma^i_{jk} = 1/2 g^{il} ( \partial_j g_{lk} + \partial_k g_{jl} - \partial_l g_{jk} )
	 * where g_{ij} is the metric tensor and g^{ij} is the inverse of the metric tensor
	 * This use the standard formula to compute the Christoffel symbols in ADM 3+1 formalism
	 * */

	for(int kk=0; kk<3; kk++){
		for(int aa=0; aa<3; aa++){
			for(int bb=0; bb<3; bb++){
				double sum=0.0;
				for(int ll=0; ll<3; ll++){

					double tmp = dgamma[aa][ll][bb] + dgamma[bb][ll][aa] - dgamma[ll][aa][bb];
					sum += invg[kk][ll] * tmp;
				}
				grid_obj.getCell(i, j, k).Christoffel[kk][aa][bb] = 0.5 * sum;
				
			}
		}
	}
	/* if (i == NX/2 && j == NY/2 && k == NZ/2) { */
	/* 	printf("Dérivées de gamma :\n"); */
	/* 	print_matrix_2D("dgamma/dx", dgamma[0]); */
	/* 	print_matrix_2D("dgamma/dy", dgamma[1]); */
	/* 	print_matrix_2D("dgamma/dz", dgamma[2]); */
	/* } */
	/* printf("Vérification de la symétrie des dérivées de g_{ij} :\n"); */
	/* for (int i = 0; i < 3; i++) { */
	/* 	for (int j = 0; j < 3; j++) { */
	/* 		printf("∂_x g[%d][%d] - ∂_y g[%d][%d] = %e\n", i, j, i, j,  */
	/* 				partialX_gamma(grid_obj, NX/2, NY/2, NZ/2, i, j) -  */
	/* 				partialY_gamma(grid_obj, NX/2, NY/2, NZ/2, i, j)); */
	/* 	} */
	/* } */
	/* if (i == NX/2 && j == NY/2 && k == NZ/2) { */
	/* 	printf("Christoffel Symboles :\n"); */
	/* 	for (int a = 0; a < 3; a++) { */
	/* 		for (int b = 0; b < 3; b++) { */
	/* 			for (int c = 0; c < 3; c++) { */
	/* 				printf("Γ[%d][%d][%d] = %e\n", a, b, c, grid_obj.getCell(i, j, k).Christoffel[a][b][c]); */
	/* 			} */
	/* 		} */
	/* 	} */
	/* } */
	/* double sum = 0.0; */
	/* int count = 0; */
	/* for (int a = 0; a < 3; a++) { */
	/* 	for (int b = 0; b < 3; b++) { */
	/* 		for (int c = 0; c < 3; c++) { */
	/* 			sum += fabs(grid_obj.getCell(i, j, k).Christoffel[a][b][c]); */
	/* 			count++; */
	/* 		} */
	/* 	} */
	/* } */
	/* printf("Moyenne absolue des Christoffel: %e\n", sum / count); */
	/*  */
	/* double max_diff = 0.0; */
	/* for (int a = 0; a < 3; a++) { */
	/* 	for (int b = 0; b < 3; b++) { */
	/* 		double diff_x = fabs(partialX_gamma(grid_obj, NX/2, NY/2, NZ/2, a, b)); */
	/* 		double diff_y = fabs(partialY_gamma(grid_obj, NX/2, NY/2, NZ/2, a, b)); */
	/* 		double diff_z = fabs(partialZ_gamma(grid_obj, NX/2, NY/2, NZ/2, a, b)); */
	/* 		max_diff = fmax(max_diff, diff_x + diff_y + diff_z); */
	/* 	} */
	/* } */
	/* printf("Max des dérivées de g_{ij}: %e\n", max_diff); */
}

