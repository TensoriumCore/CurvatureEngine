#include <Geodesics.h>

/* double compute_partial_beta(Grid &grid_obj, int i, int j, int k, int beta_index, int dim, double dx) { */
/*     double beta_p = grid_obj.getCell(i + (dim == 0), j + (dim == 1), k + (dim == 2)).beta[beta_index]; */
/*     double beta_m = grid_obj.getCell(i - (dim == 0), j - (dim == 1), k - (dim == 2)).beta[beta_index]; */
/*  */
/*     if (beta_p == 0.0 || beta_m == 0.0)  */
/*         return (grid_obj.getCell(i, j, k).beta[beta_index] - beta_m) / dx; */
/*     return (beta_p - beta_m) / (2.0 * dx); */
/* } */
/*  */
void GridTensor::compute_extrinsic_curvature(Grid &grid_obj, int i, int j, int k, \
											 double dx, double dy, double dz) {
	Grid::Cell2D &cell = grid_obj.getCell(i, j, k);

    double partialBeta[3][3];

    partialBeta[0][0] = (grid_obj.getCell(i+1, j, k).beta[0] - grid_obj.getCell(i-1, j, k).beta[0]) / (2.0 * dx);
    partialBeta[0][1] = (grid_obj.getCell(i, j+1, k).beta[0] - grid_obj.getCell(i, j-1, k).beta[0]) / (2.0 * dy);
    partialBeta[0][2] = (grid_obj.getCell(i, j, k+1).beta[0] - grid_obj.getCell(i, j, k-1).beta[0]) / (2.0 * dz);

    partialBeta[1][0] = (grid_obj.getCell(i+1, j, k).beta[1] - grid_obj.getCell(i-1, j, k).beta[1]) / (2.0 * dx);
    partialBeta[1][1] = (grid_obj.getCell(i, j+1, k).beta[1] - grid_obj.getCell(i, j-1, k).beta[1]) / (2.0 * dy);
    partialBeta[1][2] = (grid_obj.getCell(i, j, k+1).beta[1] - grid_obj.getCell(i, j, k-1).beta[1]) / (2.0 * dz);

    partialBeta[2][0] = (grid_obj.getCell(i+1, j, k).beta[2] - grid_obj.getCell(i-1, j, k).beta[2]) / (2.0 * dx);
    partialBeta[2][1] = (grid_obj.getCell(i, j+1, k).beta[2] - grid_obj.getCell(i, j-1, k).beta[2]) / (2.0 * dy);
    partialBeta[2][2] = (grid_obj.getCell(i, j, k+1).beta[2] - grid_obj.getCell(i, j, k-1).beta[2]) / (2.0 * dz);
	/* for (int ii = 0; ii < 3; ii++) { */
	/* 	for (int jj = 0; jj < 3; jj++) { */
	/* 		printf("partialBeta[%d][%d] = %e\n", ii, jj, partialBeta[ii][jj]); */
	/* 	} */
	/* } */
    compute_christoffel_3D(grid_obj, i, j, k, cell.Christoffel);
	double sumGammaBeta[3][3] = {0.0};
    for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
            for (int m = 0; m < 3; m++) {
                sumGammaBeta[ii][jj] += cell.Christoffel[ii][jj][m] * cell.beta[m];
            }
        }
    }

    double alphaLoc = cell.alpha;
    for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
            double derivPart = partialBeta[ii][jj] + partialBeta[jj][ii];
            double gammaTerm = 2.0 * sumGammaBeta[ii][jj];
            cell.K[ii][jj] = (1.0 / (2.0 * alphaLoc)) * (derivPart - gammaTerm);
        }
    }
}
