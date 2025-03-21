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
    int ip1 = i, im1 = i;
    int jp1 = j, jm1 = j;
    int kp1 = k, km1 = k;

    if (dim == 0) {
        if (i == 0) { ip1 = i+1; im1 = i; }
        else if (i == NX-1) { im1 = i-1; ip1 = i; }
        else { ip1 = i+1; im1 = i-1; }
    }
    if (dim == 1) {
        if (j == 0) { jp1 = j+1; jm1 = j; }
        else if (j == NY-1) { jm1 = j-1; jp1 = j; }
        else { jp1 = j+1; jm1 = j-1; }
    }
    if (dim == 2) {
        if (k == 0) { kp1 = k+1; km1 = k; }
        else if (k == NZ-1) { km1 = k-1; kp1 = k; }
        else { kp1 = k+1; km1 = k-1; }
    }

    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
                double Gamma_p = grid_obj.getCell(ip1, jp1, kp1).Christoffel[kk][aa][bb];
                double Gamma_m = grid_obj.getCell(im1, jm1, km1).Christoffel[kk][aa][bb];

                partialGamma[dim][kk][aa][bb] = (Gamma_p - Gamma_m) / (2.0 * d);
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
 *	R_{ab} = \partial_k \grid_obj.getCell(i, j, k).Christoffel^k_{ab} - \partial_b \grid_obj.getCell(i, j, k).Christoffel^k_{ak} + \grid_obj.getCell(i, j, k).Christoffel^l_{ab} \grid_obj.getCell(i, j, k).Christoffel^k_{lk} - \grid_obj.getCell(i, j, k).Christoffel^l_{ak} \grid_obj.getCell(i, j, k).Christoffel^k_{bl}
 *	where \grid_obj.getCell(i, j, k).Christoffel^k_{ab} is the Christoffel symbol
 *	\partial_k is the partial derivative
 *	
 *	We use the partial derivative of the Christoffel symbols to compute the Ricci tensor
 *	for the Hamiltonian constraint and the momentum constraint
 *	This tensor is used to compute the Ricci scalar and extrinsic curvature tensor in the ADM formalism
 * */


void Grid::compute_ricci_3D_conformal(Grid &grid_obj, int i, int j, int k, double Ricci[3][3]) {
    GridTensor gridTensor;
    /* gridTensor.compute_christoffel_3D(grid_obj, i, j, k, grid_obj.getCell(i, j, k).Christoffel); */

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
					term3 += grid_obj.getCell(i, j, k).Christoffel[k][a][b] * grid_obj.getCell(i, j, k).Christoffel[l][l][k]; // Correcte
					term4 += grid_obj.getCell(i, j, k).Christoffel[k][a][l] * grid_obj.getCell(i, j, k).Christoffel[l][b][k]; // Correcte
				}
			}

			Ricci[a][b] = term1 - term2 - term3 + term4;
		}
	}


}


void GridTensor::compute_ricci_conformal_factor(Grid &grid_obj, int i, int j, int k, double RicciChi[3][3]) {
    Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
    double chi = cell.chi;
    double invChi = 1.0 / chi;
    double invChi2 = invChi * invChi;

    double dChi[3];
    dChi[0] = (grid_obj.getCell(i+1,j,k).chi - grid_obj.getCell(i-1,j,k).chi) / (2.0 * DX);
    dChi[1] = (grid_obj.getCell(i,j+1,k).chi - grid_obj.getCell(i,j-1,k).chi) / (2.0 * DY);
    dChi[2] = (grid_obj.getCell(i,j,k+1).chi - grid_obj.getCell(i,j,k-1).chi) / (2.0 * DZ);

    double ddChi[3][3]; 
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            ddChi[a][b] = second_partial(grid_obj, i, j, k, a, b, [&](const Grid::Cell2D &c) {
                return c.chi;
            });
        }
    }

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double term1 = 0.0, term2 = 0.0;
            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    term1 += cell.tildgamma_inv[m][n] * ddChi[m][n]; // Δχ
                    term2 += cell.tildgamma_inv[m][n] * dChi[m] * dChi[n]; // |∇χ|^2
                }
            }

            RicciChi[a][b] =
                0.5 * invChi * (ddChi[a][b] + cell.tilde_gamma[a][b] * term1)
              - 0.25 * invChi2 * (dChi[a] * dChi[b] + cell.tilde_gamma[a][b] * term2);
        }
    }
}

void GridTensor::compute_ricci_BSSN(Grid &grid_obj, int i, int j, int k, double Ricci[3][3]) {
    double RicciTilde[3][3], RicciChi[3][3];
    grid_obj.compute_ricci_3D_conformal(grid_obj, i, j, k, RicciTilde);
    compute_ricci_conformal_factor(grid_obj, i, j, k, RicciChi);

    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            Ricci[a][b] = RicciTilde[a][b] + RicciChi[a][b];
	for (int a = 0; a < 3; a++) 
		for (int b = 0; b < 3; b++) 
			grid_obj.getCell(i, j, k).Ricci[a][b] = Ricci[a][b];
	/* #pragma omp critical */
	/* { */
	/* 	print_matrix_2D("Ricci", Ricci); */
	/* } */
}
