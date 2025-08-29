#include <Geodesics.h>

void Grid::solve_lichnerowicz(int max_iter, float tol, float dx, float dy, float dz) {
    std::vector<std::vector<std::vector<float>>> psi(NX,
        std::vector<std::vector<float>>(NY, std::vector<float>(NZ, 1.0)));

    const float inv_dx2 = 1.0/(dx*dx), inv_dy2 = 1.0/(dy*dy), inv_dz2 = 1.0/(dz*dz);
    const float factor = 1.0 / (2.0*(inv_dx2 + inv_dy2 + inv_dz2));

    for (int it = 0; it < max_iter; ++it) {
        float max_diff = 0.0;

        #pragma omp parallel for collapse(3) reduction(max:max_diff)
        for (int i = 1; i < NX - 1; ++i) {
            for (int j = 1; j < NY - 1; ++j) {
                for (int k = 1; k < NZ - 1; ++k) {
                    float lap = inv_dx2*(psi[i+1][j][k] + psi[i-1][j][k] - 2.0*psi[i][j][k])
                        + inv_dy2*(psi[i][j+1][k] + psi[i][j-1][k] - 2.0*psi[i][j][k])
                        + inv_dz2*(psi[i][j][k+1] + psi[i][j][k-1] - 2.0*psi[i][j][k]);

                    Cell2D &cell = globalGrid[i][j][k];
                    float A2 = 0.0;
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b)
                            for (int c = 0; c < 3; ++c)
                                for (int d = 0; d < 3; ++d)
                                    A2 += cell.geom.tildgamma_inv[a][c] * cell.geom.tildgamma_inv[b][d] 
                                        * cell.atilde.Atilde[c][d] * cell.atilde.Atilde[a][b];

                    float psi7 = std::pow(std::fmax(psi[i][j][k], 1e-8), 7.0);
                    float rhs = (1.0/8.0) * A2 / psi7;
                    float new_psi = psi[i][j][k] + 0.5 * (lap - rhs) * factor;
                    new_psi = std::fmax(new_psi, 1e-6f);
                    max_diff = std::fmax(max_diff, std::abs(new_psi - psi[i][j][k]));
                    psi[i][j][k] = new_psi;
                }
            }
        }

        #pragma omp parallel for collapse(3)
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int k = 0; k < NZ; k++) {
                    if (i==0 || i==NX-1 || j==0 || j==NY-1 || k==0 || k==NZ-1) {
                        float x = (i - NX/2)*dx;
                        float y = (j - NY/2)*dy;
                        float z = (k - NZ/2)*dz;
                        float r = std::sqrt(x*x+y*y+z*z);
                        psi[i][j][k] = 1.0 + 1.0/(2.0*std::max(r,1e-6f));
                    }
                }
            }
        }
		if (it % 50 == 0) {
			float minpsi = 1e9, maxpsi = -1e9;
			for (int i=0;i<NX;i++)
				for (int j=0;j<NY;j++)
					for (int k=0;k<NZ;k++) {
						minpsi = std::min(minpsi, psi[i][j][k]);
						maxpsi = std::max(maxpsi, psi[i][j][k]);
					}
			printf("Iter %d: minψ=%e maxψ=%e\n", it, minpsi, maxpsi);
		}
        if (max_diff < tol) break;
    }

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            for (int k = 0; k < NZ; ++k) {
                globalGrid[i][j][k].chi = 1.0 / std::pow(std::fmax(psi[i][j][k], 1e-8), 4);
            }
        }
    }
}


