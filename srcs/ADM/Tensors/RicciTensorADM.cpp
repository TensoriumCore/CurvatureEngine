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
    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
                auto christoffel_accessor = [&](const Grid::Cell2D &cell) {
                    return cell.Christoffel[kk][aa][bb];
                };

                if (true) { 
                    partialGamma[dim][kk][aa][bb] = second_order_diff(
                        christoffel_accessor(grid_obj.getCell(i+1, j, k)),
                        christoffel_accessor(grid_obj.getCell(i-1, j, k)),
                        d
                    );
                } else {
                    partialGamma[dim][kk][aa][bb] = fourth_order_diff(
                        christoffel_accessor(grid_obj.getCell(i+2, j, k)),
                        christoffel_accessor(grid_obj.getCell(i+1, j, k)),
                        christoffel_accessor(grid_obj.getCell(i-1, j, k)),
                        christoffel_accessor(grid_obj.getCell(i-2, j, k)),
                        d
                    );
                }
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
if (i == NX/2 && j == NY/2 && k == NZ/2) {
    printf("Christoffel Γ^0_01 = %e\n", grid_obj.getCell(i, j, k).Christoffel[0][0][1]);
    printf("Christoffel Γ^1_11 = %e\n", grid_obj.getCell(i, j, k).Christoffel[1][1][1]);
    printf("Christoffel Γ^2_12 = %e\n", grid_obj.getCell(i, j, k).Christoffel[2][1][2]);
}
    double partialGamma[3][3][3][3] = {};
	Log logger;
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 0, partialGamma, DX);
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 1, partialGamma, DY);
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 2, partialGamma, DZ);
	
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;
            for (int m = 0; m < 3; m++) {
                term1 += partialGamma[m][m][a][b];
                term2 += partialGamma[a][m][m][b];
            }

            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    term3 += Gamma[k][a][b] * Gamma[l][k][l];
                    term4 += Gamma[l][a][k] * Gamma[k][b][l];
                }
            }

            Ricci[a][b] = term1 - term2 - term3 + term4;
        }
    }

/* #pragma omp critical */
/*     { */
/*         print_matrix_2D("Ricci", Ricci); */
/*     } */
}
