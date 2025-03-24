#include <Geodesics.h>
#include <cassert>


void Grid::compute_time_derivatives(Grid &grid_obj, int i, int j, int k)
{
    Cell2D &cell = globalGrid[i][j][k];
    BSSNevolve bssn;
    double alpha = cell.gauge.alpha;
    double beta[3] = { cell.gauge.beta[0], cell.gauge.beta[1], cell.gauge.beta[2] };

      GridTensor gridTensor;
    double Gamma[3][3][3];
    gridTensor.compute_christoffel_3D(grid_obj, i, j, k, Gamma);

    double Ricci[3][3];
    gridTensor.compute_ricci_BSSN(grid_obj, i, j, k, Ricci);

    double partialBeta[3][3];
    for (int dim = 0; dim < 3; ++dim) {
        for (int comp = 0; comp < 3; ++comp) {
            partialBeta[dim][comp] = partial_m(grid_obj, i, j, k, dim,
                [&](const Grid::Cell2D &c) { return c.gauge.beta[comp]; }
            );
        }
    }

    double partialAlpha[3];
    for (int dim = 0; dim < 3; ++dim) {
        partialAlpha[dim] = partial_m(grid_obj, i, j, k, dim,
            [&](const Grid::Cell2D &c) { return c.gauge.alpha; }
        );
    }

    double d2Alpha[3][3];
    for (int m2 = 0; m2 < 3; ++m2) {
        for (int n2 = 0; n2 < 3; ++n2) {
            d2Alpha[m2][n2] = second_partial_alpha(grid_obj, i, j, k, m2, n2);
        }
    }

    double partialKtrace[3];
    for (int dim = 0; dim < 3; ++dim) {
        partialKtrace[dim] = partial_m(grid_obj, i, j, k, dim,
            [&](const Grid::Cell2D &c) { return c.curv.K_trace; }
        );
    }

    double partialTildeGamma[3][3][3];
    double partialAtilde[3][3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int dim = 0; dim < 3; dim++) {
                partialTildeGamma[dim][a][b] = partial_m(grid_obj, i, j, k, dim,
                    [&](const Grid::Cell2D &c) { return c.geom.tilde_gamma[a][b]; }
                );
                partialAtilde[dim][a][b] = partial_m(grid_obj, i, j, k, dim,
                    [&](const Grid::Cell2D &c) { return c.atilde.Atilde[a][b]; }
                );
            }
        }
    }

    double dt_chi = 0.0;
    bssn.compute_dt_chi(grid_obj, i, j, k, dt_chi);
    cell.dt_chi = dt_chi;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double adv = 0.0;
            for (int m = 0; m < 3; ++m) {
                adv += beta[m] * partialTildeGamma[m][a][b];
            }

            double shift = 0.0;
            for (int m = 0; m < 3; ++m) {
                shift += cell.geom.tilde_gamma[a][m] * partialBeta[m][b];
                shift += cell.geom.tilde_gamma[b][m] * partialBeta[m][a];
            }

            cell.geom.dt_tilde_gamma[a][b] = -2.0 * alpha * cell.atilde.Atilde[a][b] + adv + shift;
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {

            double D2_alpha = d2Alpha[a][b];
            double sumG = 0.0;
            for (int m = 0; m < 3; ++m) {
                sumG += Gamma[m][a][b] * partialAlpha[m];
            }
            D2_alpha -= sumG;

            double trace_D2_alpha = 0.0;
#pragma omp simd reduction(+:trace_D2_alpha)
            for (int mm = 0; mm < 3; ++mm) {
                for (int nn = 0; nn < 3; ++nn) {
                    double part = d2Alpha[mm][nn];
                    double gpart = 0.0;
                    for (int ll = 0; ll < 3; ++ll) {
                        gpart += Gamma[ll][mm][nn] * partialAlpha[ll];
                    }
                    trace_D2_alpha += cell.geom.tildgamma_inv[mm][nn] * (part - gpart);
                }
            }

            double Ricci_TF = Ricci[a][b];
            double R_scalar = 0.0;
            for (int mm = 0; mm < 3; ++mm) {
                for (int nn = 0; nn < 3; ++nn) {
                    R_scalar += cell.geom.tildgamma_inv[mm][nn] * Ricci[mm][nn];
                }
            }
            Ricci_TF -= (1.0/3.0) * cell.geom.tilde_gamma[a][b] * R_scalar;

            double A_A = 0.0;
            for (int k1 = 0; k1 < 3; ++k1) {
                for (int l1 = 0; l1 < 3; ++l1) {
                    A_A += cell.atilde.Atilde[a][k1]
                         * cell.geom.tildgamma_inv[k1][l1]
                         * cell.atilde.Atilde[l1][b];
                }
            }

            double adv = 0.0;
            for (int m = 0; m < 3; ++m) {
                adv += beta[m] * partialAtilde[m][a][b];
            }

            double shift_term = 0.0;
            for (int m = 0; m < 3; ++m) {
                shift_term += cell.atilde.Atilde[m][b] * partialBeta[m][a];
                shift_term += cell.atilde.Atilde[a][m] * partialBeta[m][b];
            }

            cell.atilde.dt_Atilde[a][b] =
                  cell.chi * (
                      -D2_alpha
                      + (1.0/3.0) * cell.geom.tilde_gamma[a][b] * trace_D2_alpha
                      + alpha * Ricci_TF
                  )
                + alpha * (
                      cell.curv.K_trace * cell.atilde.Atilde[a][b]
                      - 2.0 * A_A
                  )
                + adv
                + shift_term;
        }
    }
    double laplacian_alpha = 0.0;
    for (int mm = 0; mm < 3; ++mm) {
        for (int nn = 0; nn < 3; ++nn) {
            double d2a = d2Alpha[mm][nn];
            double gamma_conn = 0.0;
#pragma omp simd
            for (int ll = 0; ll < 3; ++ll) {
                gamma_conn += Gamma[ll][mm][nn] * partialAlpha[ll];
            }
            laplacian_alpha += cell.geom.tildgamma_inv[mm][nn] * (d2a - gamma_conn);
        }
    }

    double Atilde_squared = 0.0;
#pragma omp parallel for collapse(4) reduction(+:Atilde_squared)
    for (int a1 = 0; a1 < 3; ++a1) {
        for (int b1 = 0; b1 < 3; ++b1) {
            for (int c1 = 0; c1 < 3; ++c1) {
                for (int d1 = 0; d1 < 3; ++d1) {
                    Atilde_squared +=
                        cell.geom.tildgamma_inv[a1][c1] *
                        cell.geom.tildgamma_inv[b1][d1] *
                        cell.atilde.Atilde[a1][b1] *
                        cell.atilde.Atilde[c1][d1];
                }
            }
        }
    }

    double total_R_scalar = 0.0;
#pragma omp simd reduction(+:total_R_scalar)
    for (int a1 = 0; a1 < 3; ++a1) {
        for (int b1 = 0; b1 < 3; ++b1) {
            total_R_scalar += cell.geom.tildgamma_inv[a1][b1] * Ricci[a1][b1];
        }
    }

    double adv_K = 0.0;
    for (int m = 0; m < 3; ++m) {
        adv_K += beta[m] * partialKtrace[m];
    }

    cell.curv.dt_K_trace =
        -laplacian_alpha
        + alpha * (
              Atilde_squared
            + (1.0/3.0)*cell.curv.K_trace*cell.curv.K_trace
            + total_R_scalar
          )
        + adv_K;
}
