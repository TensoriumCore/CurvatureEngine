#include <Geodesics.h>

void BSSNevolve::compute_dt_tilde_gamma(Grid &grid_obj, int i, int j, int k, double dt_tg[3][3]) {
    const auto &cell = grid_obj.getCell(i, j, k);
    const double alpha = cell.gauge.alpha;

    const double* beta = cell.gauge.beta;
    const auto& tg = cell.geom.tilde_gamma;

    double d_tg[3][3][3]; 
    for (int dir = 0; dir < 3; ++dir) {
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                d_tg[dir][a][b] = partial_m(grid_obj, i, j, k, dir, [&](const Grid::Cell2D &c) {
                    return c.geom.tilde_gamma[a][b];
                });
            }
        }
    }

    double d_beta[3][3]; // j, k
    for (int j_idx = 0; j_idx < 3; ++j_idx) {
        for (int k_idx = 0; k_idx < 3; ++k_idx) {
            d_beta[j_idx][k_idx] = partial_m(grid_obj, i, j, k, j_idx, [&](const Grid::Cell2D &c) {
                return c.gauge.beta[k_idx];
            });
        }
    }

    double Lie[3][3] = {};
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            // β^k ∂_k γ̃_ij
            for (int dir = 0; dir < 3; ++dir)
                Lie[a][b] += beta[dir] * d_tg[dir][a][b];

            for (int k_idx = 0; k_idx < 3; ++k_idx)
                Lie[a][b] += tg[a][k_idx] * d_beta[b][k_idx] + tg[b][k_idx] * d_beta[a][k_idx];
        }
    }

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            dt_tg[a][b] = -2.0 * alpha * cell.atilde.Atilde[a][b] + Lie[a][b];
        }
    }
}

