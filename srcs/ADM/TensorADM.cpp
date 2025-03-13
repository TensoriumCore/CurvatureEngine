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


void GridTensor::compute_partial_christoffel(Grid &grid_obj, int i, int j, int k, int dim, double partialGamma[3][3][3][3], double d) {
    double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];
    double localPartialGamma[3][3][3] = {0.0}; 

    int a = (dim == 0) ? i : ((dim == 1) ? j : k);
    int max_a = (dim == 0) ? NX : ((dim == 1) ? NY : NZ);

    if (a >= 2 && a <= max_a - 3) {
        compute_christoffel_3D(grid_obj, i - 2 * (dim == 0), j - 2 * (dim == 1), k - 2 * (dim == 2), Gmm);
        compute_christoffel_3D(grid_obj, i - 1 * (dim == 0), j - 1 * (dim == 1), k - 1 * (dim == 2), Gm);
        compute_christoffel_3D(grid_obj, i + 1 * (dim == 0), j + 1 * (dim == 1), k + 1 * (dim == 2), Gp);
        compute_christoffel_3D(grid_obj, i + 2 * (dim == 0), j + 2 * (dim == 1), k + 2 * (dim == 2), Gpp);
    } else if (a >= 1 && a <= max_a - 2) {
        compute_christoffel_3D(grid_obj, i - 1 * (dim == 0), j - 1 * (dim == 1), k - 1 * (dim == 2), Gm);
        compute_christoffel_3D(grid_obj, i + 1 * (dim == 0), j + 1 * (dim == 1), k + 1 * (dim == 2), Gp);
    }

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
                partialGamma[dim][kk][aa][bb] = localPartialGamma[kk][aa][bb];
            }
        }
    }
}


void Grid::compute_ricci_3D(Grid &grid_obj, int i, int j, int k, double Ricci[3][3]) {
    double Gamma[3][3][3];
    GridTensor gridTensor;
    gridTensor.compute_christoffel_3D(grid_obj, i, j, k, Gamma);

    double T[3][3]; 
    double partialGamma[3][3][3][3] = {};  
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 0, partialGamma, DX);
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 1, partialGamma, DY);
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 2, partialGamma, DZ);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;
            
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

            Ricci[a][b] = term1 - term2 - term3 + term4;
        }
    }
	/* if (i == NX/2 && j == NY/2 && k == NZ/2) { */
	/* 	print_matrix_2D("Gamma", Gamma[0]); */
	/* } */
	/* print_matrix_2D("Ricci", Ricci); */
}
