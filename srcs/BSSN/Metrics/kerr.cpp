#include <Geodesics.h>

static inline float kerrSchildRadius(float x, float y, float z, float a)
{
    const float r2   = x*x + y*y + z*z;
    const float a2   = a*a;
    const float alpha = r2 - a2;
    const float inside = alpha*alpha + 4.0*a2*z*z;

    const float term = 0.5 * (alpha + std::sqrt(inside));

    if (term <= 0.0) {
        return 0.0;
    }
    return std::sqrt(term);
}

void Grid::initializeKerrData(Grid &grid_obj) {
    float m = 1.0;
    float a = 0.9999;
    float x0 = 0.0, y0 = 0.0, z0 = 0.0;

    float L = 9.0;
    float x_min = -L, x_max = L;
    float y_min = -L, y_max = L;
    float z_min = -L, z_max = L;

    float dx = (x_max - x_min) / (NX - 1);
    float dy = (y_max - y_min) / (NY - 1);
    float dz = (z_max - z_min) / (NZ - 1);

    GridTensor gridtensor;
    Matrix matrix;

    globalGrid.resize(NX, std::vector<std::vector<Cell2D>>(NY, std::vector<Cell2D>(NZ)));

    #pragma omp parallel for collapse(3) schedule(dynamic)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                float x = x_min + i * dx;
                float y = y_min + j * dy;
                float z = z_min + k * dz;

                float dx0 = x - x0;
                float dy0 = y - y0;
                float dz0 = z - z0;

				float rKS = kerrSchildRadius(dx0, dy0, dz0, a);
                float cosTheta = (rKS > 1e-14) ? dz0 / rKS : 0.0;
                float denomKS = (rKS * rKS) + (a * a * cosTheta * cosTheta);
                float H = (rKS > 1e-14 && denomKS > 1e-14) ? (m * rKS) / denomKS : 0.0;

                float denomVec = (rKS * rKS) + (a * a);
                float lx = (denomVec > 1e-14) ? (rKS * dx0 + a * dy0) / denomVec : 0.0;
                float ly = (denomVec > 1e-14) ? (rKS * dy0 - a * dx0) / denomVec : 0.0;
                float lz = (rKS > 1e-14) ? dz0 / rKS : 0.0;

                float norm_l = std::sqrt(lx * lx + ly * ly + lz * lz);
                if (norm_l > 1e-14) {
                    lx /= norm_l;
                    ly /= norm_l;
                    lz /= norm_l;
                }

                Cell2D &cell = globalGrid[i][j][k];

                cell.geom.gamma[0][0] = 1.0 + 2.0 * H * lx * lx;
                cell.geom.gamma[0][1] = 2.0 * H * lx * ly;
                cell.geom.gamma[0][2] = 2.0 * H * lx * lz;
                cell.geom.gamma[1][1] = 1.0 + 2.0 * H * ly * ly;
                cell.geom.gamma[1][2] = 2.0 * H * ly * lz;
                cell.geom.gamma[2][2] = 1.0 + 2.0 * H * lz * lz;
                cell.geom.gamma[1][0] = cell.geom.gamma[0][1];
                cell.geom.gamma[2][0] = cell.geom.gamma[0][2];
                cell.geom.gamma[2][1] = cell.geom.gamma[1][2];

                matrix.inverse_3x3(cell.geom.gamma, cell.geom.gamma_inv);

                float det_gamma = std::cbrt(
                    cell.geom.gamma[0][0] * cell.geom.gamma[1][1] * cell.geom.gamma[2][2]
                  - cell.geom.gamma[0][0] * cell.geom.gamma[1][2] * cell.geom.gamma[2][1]
                  - cell.geom.gamma[1][1] * cell.geom.gamma[0][2] * cell.geom.gamma[2][0]
                  - cell.geom.gamma[2][2] * cell.geom.gamma[0][1] * cell.geom.gamma[1][0]
                );
                cell.chi = 1.0 / det_gamma;

                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        cell.geom.tilde_gamma[a][b] = cell.chi * cell.geom.gamma[a][b];

                Matrix3x3 tilde_gamma_std;
                Matrix3x3 tilde_gamma_inv_std;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        tilde_gamma_std[a][b] = cell.geom.tilde_gamma[a][b];
                matrix.inverse_3x3(tilde_gamma_std, tilde_gamma_inv_std);
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        cell.geom.tildgamma_inv[a][b] = tilde_gamma_inv_std[a][b];

                cell.gauge.alpha = 1.0 / std::sqrt(1.0 + 2.0 * H);
                cell.gauge.beta[0] = 2.0 * H * lx;
                cell.gauge.beta[1] = 2.0 * H * ly;
                cell.gauge.beta[2] = 2.0 * H * lz;

                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        cell.curv.K[a][b] = 0.0;
            }
        }
    }

    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int k = 1; k < NZ - 1; k++) {
                gridtensor.compute_extrinsic_curvature(*this, i, j, k, dx, dy, dz);
                gridtensor.compute_Atilde(*this, i, j, k);
            }
        }
    }
#ifdef FLUID
	initializeFishboneMoncriefTorus(3.0, 6.0, 1.0, 3.8, 4.0);
#endif
	printf("K_{ij} near (1,1,1) = \n");
    for (int p = 0; p < 3; p++) {
        printf("  ");
        for (int q = 0; q < 3; q++) {
            printf("%e ", globalGrid[1][1][1].curv.K[p][q]);
        }
        printf("\n");
    }
    printf("Gamma_{ij} near (1,1,1) = \n");
    for (int p = 0; p < 3; p++) {
        printf("  ");
        for (int q = 0; q < 3; q++) {
            printf("%e ", globalGrid[1][1][1].geom.gamma[p][q]);
        }
        printf("\n");
    }
    printf("chi near (1,1,1) = %e\n", globalGrid[1][1][1].chi);
    printf("alpha_1_1_1 = %e\n", globalGrid[1][1][1].gauge.alpha);
    
    float test_radii[] = {0.5, 1.0, 2.0, 5.0, 10.0};
    for (float test_r : test_radii) {
        float H_test = M / test_r;  
        float alpha_test = 1.0 / std::sqrt(1.0 + 2.0 * H_test);
        printf("Eq plane r = %f : H = %e, alpha = %e (Schw approx)\n", test_r, H_test, alpha_test);
    }

}

