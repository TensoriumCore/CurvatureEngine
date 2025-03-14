#include <Geodesics.h>

static void print_matrix_2D(const char* name, const double mat[3][3])
{
    printf("\n%s (3x3):\n", name);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%12.6f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/*
 * This function compute the partial derivative of the Christoffel symbols at each point of the 3D grid
 * @param i , j , k the index of the cell
 * @param dim the dimension of the partial derivative
 * @param partialGamma the output array of the partial derivative of the Christoffel tensor
 * @param d the step size
 * @return void
 *
 * this function is based on the fourth order finite difference method
 * it currently only works for the 3D case and need to be analyticaly tested
 * */

void GridTensor::compute_partial_christoffel(Grid &grid_obj, int i, int j, int k, int dim, double partialGamma[3][3][3][3], double d) {
    double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];
    double localPartialGamma[3][3][3] = {0.0}; 

    int a = (dim == 0) ? i : ((dim == 1) ? j : k);
    int max_a = (dim == 0) ? NX : ((dim == 1) ? NY : NZ);

	/*
	 * This is the fourth order finite difference method
	 * we compute the partial derivative of the Christoffel symbols using the fourth order finite difference method
	 * */

    if (a >= 2 && a <= max_a - 3) {
        compute_christoffel_3D(grid_obj, i - 2 * (dim == 0), j - 2 * (dim == 1), k - 2 * (dim == 2), Gmm);
        compute_christoffel_3D(grid_obj, i - 1 * (dim == 0), j - 1 * (dim == 1), k - 1 * (dim == 2), Gm);
        compute_christoffel_3D(grid_obj, i + 1 * (dim == 0), j + 1 * (dim == 1), k + 1 * (dim == 2), Gp);
        compute_christoffel_3D(grid_obj, i + 2 * (dim == 0), j + 2 * (dim == 1), k + 2 * (dim == 2), Gpp);
    } else if (a >= 1 && a <= max_a - 2) {
        compute_christoffel_3D(grid_obj, i - 1 * (dim == 0), j - 1 * (dim == 1), k - 1 * (dim == 2), Gm);
        compute_christoffel_3D(grid_obj, i + 1 * (dim == 0), j + 1 * (dim == 1), k + 1 * (dim == 2), Gp);
    }

	/*
	 * This is the fourth order finite difference method
	 * If the point is not on the boundary we use the fourth order finite difference method
	 * If the point is on the boundary we use the second order finite difference method
	 * else we set the partial derivative to zero
	 * Not really sure if this is the best way to handle the boundary about the partial derivative of the Christoffel symbols
	 * */

    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
                if (a >= 2 && a <= max_a - 3) {
                    localPartialGamma[kk][aa][bb] = fourth_order_diff(Gpp[kk][aa][bb],  \
																		Gp[kk][aa][bb], \
																		Gm[kk][aa][bb], \
																		Gmm[kk][aa][bb], d);
                } else if (a >= 1 && a <= max_a - 2) {
                    localPartialGamma[kk][aa][bb] = second_order_diff(Gp[kk][aa][bb], \
																		Gm[kk][aa][bb], d);
                } else {
                    localPartialGamma[kk][aa][bb] = 0.0;
                }
				/* print_matrix_2D("Partialgamma", localPartialGamma[0]); */
            }
        }
    }

    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
				
				/*
				 * We set the partial derivative of the Christoffel symbols in the output array
				 * the formula is:
				 * \partial_a \Gamma^k_{ab} = \partial_a \Gamma^k_{ab}
				 * */

                partialGamma[dim][kk][aa][bb] = localPartialGamma[kk][aa][bb];
            }
        }
    }
}

/*
 * This function is used to compute the 3D ADM Ricci tensor
 *	@params grid_obj the grid object
 *	@params i, j, k the index of the cell
 *	@params Ricci the output array of the Ricci tensor
 *	@return void
 *
 *	The Ricci tensor is computed using the formula:
 *	R_{ab} = \partial_k \Gamma^k_{ab} - \partial_b \Gamma^k_{ak} + \Gamma^l_{ab} \Gamma^k_{lk} - \Gamma^l_{ak} \Gamma^k_{bl}
 *	where \Gamma^k_{ab} is the Christoffel symbol
 *	\partial_k is the partial derivative
 *	
 *	We use the partial derivative of the Christoffel symbols to compute the Ricci tensor
 *	for the Hamiltonian constraint and the momentum constraint
 *	This tensor is used to compute the Ricci scalar and extrinsic curvature tensor in the ADM formalism
 * */


void Grid::compute_ricci_3D(Grid &grid_obj, int i, int j, int k, double Ricci[3][3]) {
    double Gamma[3][3][3];
    GridTensor gridTensor;
    gridTensor.compute_christoffel_3D(grid_obj, i, j, k, Gamma);

    double T[3][3]; 
    double partialGamma[3][3][3][3] = {}; 

	/*
	 * We compute the partial derivative of the Christoffel symbols
	 * using the fourth order finite difference method
	 * */

    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 0, partialGamma, DX);
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 1, partialGamma, DY);
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 2, partialGamma, DZ);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;
            /*
			 * We compute the Ricci tensor using the formula:
			 * R_{ab} = \partial_k \Gamma^k_{ab} - \partial_b \Gamma^k_{ak} + \Gamma^l_{ab} \Gamma^k_{lk} - \Gamma^l_{ak} \Gamma^k_{bl}
			 * the terms are computed using the Christoffel symbols
			 * then we set the Ricci tensor in the output array
			 * */
            for (int m = 0; m < 3; m++) {
                term1 += partialGamma[m][m][a][b];  
                term2 += partialGamma[a][m][m][b];
				/* printf("term1: %f, term2: %f\n", term1, term2); */
            }
			#pragma omp simd collapse(2)
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    term3 += Gamma[k][a][b] * Gamma[l][k][l]; // Γ^k_ab Γ^l_kl
                    term4 += Gamma[l][a][k] * Gamma[k][b][l]; // Γ^l_ak Γ^k_bl
                }
            }
			
			/*
			 * We set the Ricci tensor in the output array
			 * by assigning the value of the terms in the formula
			 * */

            Ricci[a][b] = term1 - term2 - term3 + term4;
        }
    }
	/* if (i == NX/2 && j == NY/2 && k == NZ/2) { */
	/* 	print_matrix_2D("Gamma", Gamma[0]); */
	/* } */
	/* print_matrix_2D("Ricci", Ricci); */
}
