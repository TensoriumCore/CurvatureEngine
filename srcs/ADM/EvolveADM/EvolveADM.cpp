#include <Geodesics.h>
#include <cassert>


void Grid::compute_time_derivatives(Grid &grid_obj, int i, int j, int k)
{
    Cell2D &cell = globalGrid[i][j][k];
    BSSNevolve bssn;
    double alpha = cell.alpha;
    double beta[3] = { cell.beta[0], cell.beta[1], cell.beta[2] };

    GridTensor gridTensor;
    double Gamma[3][3][3];
    gridTensor.compute_christoffel_3D(grid_obj, i, j, k, Gamma);

    double Ricci[3][3];
    gridTensor.compute_ricci_BSSN(grid_obj, i, j, k, Ricci);

    double partialBeta[3][3];
    for (int jDim = 0; jDim < 3; ++jDim)
        for (int iComp = 0; iComp < 3; ++iComp)
            partialBeta[jDim][iComp] = partial_m(grid_obj, i, j, k, jDim, [&](const Grid::Cell2D &c) {
                return c.beta[iComp];
            });

    double dt_chi = 0.0;
    bssn.compute_dt_chi(grid_obj, i, j, k, dt_chi);
    cell.dt_chi = dt_chi;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double adv = 0.0;
            for (int m = 0; m < 3; ++m)
                adv += beta[m] * partial_m(grid_obj, i, j, k, m, [&](const Grid::Cell2D &c) {
                    return c.tilde_gamma[a][b];
                });

            double shift = 0.0;
            for (int m = 0; m < 3; ++m) {
                shift += cell.tilde_gamma[a][m] * partialBeta[b][m];
                shift += cell.tilde_gamma[b][m] * partialBeta[a][m];
            }

            cell.dt_tilde_gamma[a][b] = -2.0 * alpha * cell.Atilde[a][b] + adv + shift;
        }
    }

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double D2_alpha = second_partial_alpha(grid_obj, i, j, k, a, b);
            double sumG = 0.0;
            for (int m = 0; m < 3; ++m)
                sumG += Gamma[m][a][b] * partial_m(grid_obj, i, j, k, m, [&](const Grid::Cell2D &c) { return c.alpha; });
            D2_alpha -= sumG;

            double trace_D2_alpha = 0.0;
            for (int m = 0; m < 3; ++m)
                for (int n = 0; n < 3; ++n) {
                    double part = second_partial_alpha(grid_obj, i, j, k, m, n);
                    double gpart = 0.0;
                    for (int l = 0; l < 3; ++l)
                        gpart += Gamma[l][m][n] * partial_m(grid_obj, i, j, k, l, [&](const Grid::Cell2D &c) { return c.alpha; });
                    trace_D2_alpha += cell.tildgamma_inv[m][n] * (part - gpart);
                }

            double Ricci_TF = Ricci[a][b];
            double R_scalar = 0.0;
            for (int m = 0; m < 3; ++m)
                for (int n = 0; n < 3; ++n)
                    R_scalar += cell.tildgamma_inv[m][n] * Ricci[m][n];
            Ricci_TF -= (1.0 / 3.0) * cell.tilde_gamma[a][b] * R_scalar;

            double A_A = 0.0;
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    A_A += cell.Atilde[a][k] * cell.tildgamma_inv[k][l] * cell.Atilde[l][b];

            double adv = 0.0;
            for (int m = 0; m < 3; ++m)
                adv += beta[m] * partial_m(grid_obj, i, j, k, m, [&](const Grid::Cell2D &c) {
                    return c.Atilde[a][b];
                });

            double shift_term = 0.0;
            for (int m = 0; m < 3; ++m) {
                shift_term += cell.Atilde[m][b] * partialBeta[a][m];
                shift_term += cell.Atilde[a][m] * partialBeta[b][m];
            }

            cell.dt_Atilde[a][b] = cell.chi * (-D2_alpha + (1.0 / 3.0) * cell.tilde_gamma[a][b] * trace_D2_alpha + alpha * Ricci_TF)
                                 + alpha * (cell.K_trace * cell.Atilde[a][b] - 2.0 * A_A)
                                 + adv + shift_term;

        }
    }

    double laplacian_alpha = 0.0;
    for (int m = 0; m < 3; ++m) {
        for (int n = 0; n < 3; ++n) {
            double d2_alpha = second_partial_alpha(grid_obj, i, j, k, m, n);
            double gamma_conn = 0.0;
            for (int l = 0; l < 3; ++l)
                gamma_conn += Gamma[l][m][n] * partial_m(grid_obj, i, j, k, l, [&](const Grid::Cell2D &c) { return c.alpha; });

            laplacian_alpha += cell.tildgamma_inv[m][n] * (d2_alpha - gamma_conn);
        }
    }

    double Atilde_squared = 0.0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            for (int c = 0; c < 3; ++c)
                for (int d = 0; d < 3; ++d)
                    Atilde_squared += cell.tildgamma_inv[a][c] * cell.tildgamma_inv[b][d] * cell.Atilde[a][b] * cell.Atilde[c][d];

    double R_scalar = 0.0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            R_scalar += cell.tildgamma_inv[a][b] * Ricci[a][b];

    double adv_K = 0.0;
    for (int m = 0; m < 3; ++m)
        adv_K += beta[m] * partial_m(grid_obj, i, j, k, m, [&](const Grid::Cell2D &c) { return c.K_trace; });

    cell.dt_K_trace = -laplacian_alpha + alpha * (Atilde_squared + (1.0 / 3.0) * cell.K_trace * cell.K_trace + R_scalar) + adv_K;
}
