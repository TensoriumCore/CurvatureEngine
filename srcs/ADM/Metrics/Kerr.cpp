#include <Geodesics.h>


inline double kerrSchildRadius(double x, double y, double z, double a)
{
    const double r2   = x*x + y*y + z*z;
    const double a2   = a*a;
    const double alpha = r2 - a2;
    const double inside = alpha*alpha + 4.0*a2*z*z;

    const double term = 0.5 * (alpha + std::sqrt(inside));

    if (term <= 0.0) {
        return 0.0;
    }
    return std::sqrt(term);
}

void Grid::initializeKerrData(Grid &grid_obj) {
    double a = 0.9999;
    double L = 9.0; 
    double x_min = -L, x_max =  L;
    double y_min = -L, y_max =  L;
    double z_min = -L, z_max =  L;

    double dx = (x_max - x_min) / (NX - 1);
    double dy = (y_max - y_min) / (NY - 1);
    double dz = (z_max - z_min) / (NZ - 1);

    GridTensor gridtensor;
    Matrix matrix;   
    std::vector<std::array<double, 3>> horizonPoints;

    globalGrid.resize(NX, std::vector<std::vector<Cell2D>>(NY, std::vector<Cell2D>(NZ)));

#pragma omp parallel for collapse(3) schedule(dynamic)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                double z = z_min + k * dz;

                double rKS = kerrSchildRadius(x, y, z, a);
                double cosTheta = (rKS > 1e-14) ? z / rKS : 0.0;
                double denomKS = (rKS * rKS) + (a * a * cosTheta * cosTheta);
                double H = (rKS > 1e-14 && denomKS > 1e-14) ? (M * rKS) / denomKS : 0.0;

                double denomVec = (rKS * rKS) + (a * a);
                double lx = (denomVec > 1e-14) ? (rKS * x + a * y) / denomVec : 0.0;
                double ly = (denomVec > 1e-14) ? (rKS * y - a * x) / denomVec : 0.0;
                double lz = (rKS > 1e-14) ? z / rKS : 0.0;

                Cell2D &cell = globalGrid[i][j][k];

                cell.gamma[0][0] = 1.0 + 2.0 * H * lx * lx;
                cell.gamma[0][1] = 2.0 * H * lx * ly;
                cell.gamma[0][2] = 2.0 * H * lx * lz;
                cell.gamma[1][1] = 1.0 + 2.0 * H * ly * ly;
                cell.gamma[1][2] = 2.0 * H * ly * lz;
                cell.gamma[2][2] = 1.0 + 2.0 * H * lz * lz;
                cell.gamma[1][0] = cell.gamma[0][1];
                cell.gamma[2][0] = cell.gamma[0][2];
                cell.gamma[2][1] = cell.gamma[1][2];

                matrix.inverse_3x3(cell.gamma, cell.gamma_inv);

				double det_gamma = std::cbrt(
						cell.gamma[0][0] * cell.gamma[1][1] * cell.gamma[2][2]
						- cell.gamma[0][0] * cell.gamma[1][2] * cell.gamma[2][1]
						- cell.gamma[1][1] * cell.gamma[0][2] * cell.gamma[2][0]
						- cell.gamma[2][2] * cell.gamma[0][1] * cell.gamma[1][0]
						);
				cell.chi = 1.0 / det_gamma;

				for (int a = 0; a < 3; a++) {
					for (int b = 0; b < 3; b++) {
						cell.tilde_gamma[a][b] = cell.chi * cell.gamma[a][b];
					}
				}

				Matrix3x3 tilde_gamma_std;
				Matrix3x3 tilde_gamma_inv_std;

				for (int a = 0; a < 3; ++a)
					for (int b = 0; b < 3; ++b)
						tilde_gamma_std[a][b] = cell.tilde_gamma[a][b];

				matrix.inverse_3x3(tilde_gamma_std, tilde_gamma_inv_std);

				for (int a = 0; a < 3; ++a)
					for (int b = 0; b < 3; ++b)
						cell.tildgamma_inv[a][b] = tilde_gamma_inv_std[a][b];

               
                cell.alpha = 1.0 / std::sqrt(1.0 + 2.0 * H);
                cell.beta[0] = 2.0 * H * lx;
                cell.beta[1] = 2.0 * H * ly;
                cell.beta[2] = 2.0 * H * lz;

                for (int a_idx = 0; a_idx < 3; a_idx++) {
                    for (int b_idx = 0; b_idx < 3; b_idx++) {
                        cell.K[a_idx][b_idx] = 0.0;
                    }
                }
            }
        }
    }

    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int k = 1; k < NZ - 1; k++) {
                gridtensor.compute_extrinsic_curvature(*this, i, j, k, dx, dy, dz);
				gridtensor.compute_Atilde(*this, i, j, k);
                Cell2D &cell = globalGrid[i][j][k];
				/* if (i == NX/2 && j == NY/2 && k == NZ/2) { */
				/* 	printf("Atilde_{ij} =\n"); */
				/* 	for (int a = 0; a < 3; ++a) { */
				/* 		for (int b = 0; b < 3; ++b) { */
				/* 			printf("  %.5e", cell.Atilde[a][b]); */
				/* 		} */
				/* 		printf("\n"); */
				/* 	} */
				/* } */
				double trace_Atilde = 0.0;
				for (int a = 0; a < 3; ++a) {
					for (int b = 0; b < 3; ++b) {
						trace_Atilde += globalGrid[1][1][1].tildgamma_inv[a][b] * globalGrid[1][1][1].Atilde[a][b];
					}
				}
                double Ktrace = 0.0;
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        Ktrace += globalGrid[i][j][k].gamma_inv[a][b] * globalGrid[i][j][k].K[a][b];
                    }
                }
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        globalGrid[i][j][k].Atilde[a][b] = globalGrid[i][j][k].K[a][b] 
                            - (1.0 / 3.0) * Ktrace * globalGrid[i][j][k].gamma[a][b];
                    }
                }
            }
        }
	}
	printf("K_{ij} near (1,1,1) = \n");
	for (int p = 0; p < 3; p++) {
		printf("  ");
		for (int q = 0; q < 3; q++) {
			printf("%e ", globalGrid[1][1][1].K[p][q]);
		}
		printf("\n");
	}
	printf("Gamma_{ij} near (1,1,1) = \n");
	for (int p = 0; p < 3; p++) {
		printf("  ");
		for (int q = 0; q < 3; q++) {
			printf("%e ", globalGrid[1][1][1].gamma[p][q]);
		}
		printf("\n");
	}
	printf("Gamma_inv_{ij} near (1,1,1) = \n");
	for (int p = 0; p < 3; p++) {
		printf("  ");
		for (int q = 0; q < 3; q++) {
			printf("%e ", globalGrid[1][1][1].gamma_inv[p][q]);
		}
		printf("\n");
	}
	printf("Gamma_tile_{ij} near (1,1,1) = \n");
	for (int p = 0; p < 3; p++) {
		printf("  ");
		for (int q = 0; q < 3; q++) {
			printf("%e ", globalGrid[1][1][1].tilde_gamma[p][q]);
		}
		printf("\n");
	}
	printf("Gamma_tilde_inv_{ij} near (1,1,1) = \n");
	for (int p = 0; p < 3; p++) {
		printf("  ");
		for (int q = 0; q < 3; q++) {
			printf("%e ", globalGrid[1][1][1].tildgamma_inv[p][q]);
		}
		printf("\n");
	}

	printf("chi near (1,1,1) = %e\n", globalGrid[1][1][1].chi);

	printf("alpha_1_1_1 = %e\n", globalGrid[1][1][1].alpha);

	double test_radii[] = {0.5, 1.0, 2.0, 5.0, 10.0};
	for (double test_r : test_radii) {
		double H_test = M / test_r;  
		double alpha_test = 1.0 / std::sqrt(1.0 + 2.0 * H_test);
		printf("Eq plane r = %f : H = %e, alpha = %e (Schw approx)\n", 
				test_r, H_test, alpha_test);
	}
}
