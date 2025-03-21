#include <Geodesics.h>

void BSSNevolve::compute_dt_chi(Grid &grid_obj, int i, int j, int k, double &dt_chi) {
    const auto &cell = grid_obj.getCell(i, j, k);

    double chi = cell.chi;
    double alpha = cell.alpha;
    double Ktrace = 0.0;

    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            Ktrace += cell.tildgamma_inv[a][b] * cell.K[a][b];

    double div_beta = 0.0;
    for (int a = 0; a < 3; ++a) {
        div_beta += partial_m(grid_obj, i, j, k, a, [&](const Grid::Cell2D &c) {
            return c.beta[a];
        });
    }

    double beta_grad_chi = 0.0;
    for (int a = 0; a < 3; ++a) {
        double d_chi = partial_m(grid_obj, i, j, k, a, [&](const Grid::Cell2D &c) {
            return c.chi;
        });
        beta_grad_chi += cell.beta[a] * d_chi;
    }

    dt_chi = (2.0 / 3.0) * chi * (alpha * Ktrace - div_beta) + beta_grad_chi;
}

