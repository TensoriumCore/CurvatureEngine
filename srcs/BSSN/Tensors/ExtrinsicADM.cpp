#include <Geodesics.h>

/* float compute_partial_beta(Grid &grid_obj, int i, int j, int k, int beta_index, int dim, float dx) { */
/*     float beta_p = grid_obj.getCell(i + (dim == 0), j + (dim == 1), k + (dim == 2)).gauge.beta[beta_index]; */
/*     float beta_m = grid_obj.getCell(i - (dim == 0), j - (dim == 1), k - (dim == 2)).gauge.beta[beta_index]; */
/*  */
/*     if (beta_p == 0.0 || beta_m == 0.0)  */
/*         return (grid_obj.getCell(i, j, k).gauge.beta[beta_index] - beta_m) / dx; */
/*     return (beta_p - beta_m) / (2.0 * dx); */
/* } */
/*  */
void GridTensor::compute_extrinsic_curvature(Grid &grid_obj, int i, int j, int k,
                                             float dx, float dy, float dz) {
    Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
    float partialBeta[3][3];
    partialBeta[0][0] = (grid_obj.getCell(i+1, j, k).gauge.beta[0] - grid_obj.getCell(i-1, j, k).gauge.beta[0]) / (2.0 * dx);
    partialBeta[0][1] = (grid_obj.getCell(i, j+1, k).gauge.beta[0] - grid_obj.getCell(i, j-1, k).gauge.beta[0]) / (2.0 * dy);
    partialBeta[0][2] = (grid_obj.getCell(i, j, k+1).gauge.beta[0] - grid_obj.getCell(i, j, k-1).gauge.beta[0]) / (2.0 * dz);

    partialBeta[1][0] = (grid_obj.getCell(i+1, j, k).gauge.beta[1] - grid_obj.getCell(i-1, j, k).gauge.beta[1]) / (2.0 * dx);
    partialBeta[1][1] = (grid_obj.getCell(i, j+1, k).gauge.beta[1] - grid_obj.getCell(i, j-1, k).gauge.beta[1]) / (2.0 * dy);
    partialBeta[1][2] = (grid_obj.getCell(i, j, k+1).gauge.beta[1] - grid_obj.getCell(i, j, k-1).gauge.beta[1]) / (2.0 * dz);

    partialBeta[2][0] = (grid_obj.getCell(i+1, j, k).gauge.beta[2] - grid_obj.getCell(i-1, j, k).gauge.beta[2]) / (2.0 * dx);
    partialBeta[2][1] = (grid_obj.getCell(i, j+1, k).gauge.beta[2] - grid_obj.getCell(i, j-1, k).gauge.beta[2]) / (2.0 * dy);
    partialBeta[2][2] = (grid_obj.getCell(i, j, k+1).gauge.beta[2] - grid_obj.getCell(i, j, k-1).gauge.beta[2]) / (2.0 * dz);

    compute_christoffel_3D(grid_obj, i, j, k, cell.conn.Christoffel);

    float GammaBeta[3][3] = {0.0};
    for (int i_idx = 0; i_idx < 3; ++i_idx) {
        for (int j_idx = 0; j_idx < 3; ++j_idx) {
            for (int k_idx = 0; k_idx < 3; ++k_idx) {
                GammaBeta[i_idx][j_idx] += cell.conn.Christoffel[i_idx][j_idx][k_idx] * cell.gauge.beta[k_idx];
            }
        }
    }

    float alpha = cell.gauge.alpha;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            float sym_grad_beta = partialBeta[a][b] + partialBeta[b][a];
            float correction = sym_grad_beta - 2.0 * GammaBeta[a][b];
            cell.curv.K[a][b] = -0.5 / alpha * (cell.dgt[a][b] - correction);
        }
    }
}
