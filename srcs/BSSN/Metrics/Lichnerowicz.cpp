#include <Geodesics.h>


void Grid::solve_lichnerowicz(int max_iter, double tol, double dx, double dy, double dz) {
    std::vector<std::vector<std::vector<double>>> psi(NX,
        std::vector<std::vector<double>>(NY, std::vector<double>(NZ, 1.0)));

    const double inv_dx2 = 1.0/(dx*dx), inv_dy2 = 1.0/(dy*dy), inv_dz2 = 1.0/(dz*dz);
    const double factor = 1.0 / (2.0*(inv_dx2 + inv_dy2 + inv_dz2));

    for (int it = 0; it < max_iter; ++it) {
        double max_diff = 0.0;


        #pragma omp parallel for collapse(3) reduction(max:max_diff)
        for (int i = 1; i < NX - 1; ++i) {
            for (int j = 1; j < NY - 1; ++j) {
                for (int k = 1; k < NZ - 1; ++k) {
                    double lap = inv_dx2*(psi[i+1][j][k] + psi[i-1][j][k] - 2.0*psi[i][j][k])
                                + inv_dy2*(psi[i][j+1][k] + psi[i][j-1][k] - 2.0*psi[i][j][k])
                                + inv_dz2*(psi[i][j][k+1] + psi[i][j][k-1] - 2.0*psi[i][j][k]);

                    Cell2D &cell = globalGrid[i][j][k];
                    double A2 = 0.0;
                    for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            double temp = 0.0;
                            for (int c = 0; c < 3; ++c) {
                                temp += cell.geom.tildgamma_inv[a][c] * cell.atilde.Atilde[c][b];
                            }
                            A2 += temp * cell.atilde.Atilde[a][b];
                        }
                    }

                    double psi7 = std::pow(std::max(psi[i][j][k], 1e-8), 7.0);
                    double rhs = (1.0/8.0) * A2 / psi7;

                    double new_psi = psi[i][j][k] + 0.5 * (lap - rhs) * factor;
                    new_psi = std::max(new_psi, 0.1); // Garantit ψ > 0

                    max_diff = std::max(max_diff, std::abs(new_psi - psi[i][j][k]));
                    psi[i][j][k] = new_psi;
                }
            }
        }

        if (max_diff < tol) break;
    }

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            for (int k = 0; k < NZ; ++k) {
                globalGrid[i][j][k].chi = 1.0 / std::pow(std::max(psi[i][j][k], 1e-8), 4);
            }
        }
    }
}
