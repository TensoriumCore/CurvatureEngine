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

void Grid::injectTTWave(Cell2D &cell, double x, double y, double z, double t){
	constexpr double A      = 1.6; 
	constexpr double lambda = 3.0; 
	constexpr double r0     = 2.0;
	constexpr double sigma  = 0.6; 
	constexpr double omega  = 2.0 * M_PI / lambda;
	
	double r = std::sqrt(x * x + y * y + z * z);
	if (r < 1e-14) return; 

	double theta = std::acos(z / r);
	double phi   = std::atan2(y, x);

	double envelope = std::exp(-((r - r0) * (r - r0)) / (sigma * sigma));
	double phase    = omega * t - 2.0 * phi; 

	double h_plus   = A / r * std::cos(phase) * envelope;
	double h_cross  = A / r * std::sin(phase) * envelope;

	double sin_theta = std::sin(theta), cos_theta = std::cos(theta);
	double sin_phi   = std::sin(phi),   cos_phi   = std::cos(phi);

	Vector3 e_theta = {
		cos_theta * cos_phi,
		cos_theta * sin_phi,
		-sin_theta
	};

	Vector3 e_phi = {
		-sin_phi,
		cos_phi,
		0.0
	};

	Matrix3x3 h_cartesian = {0};

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) {
			h_cartesian[i][j] += h_plus * (e_theta[i] * e_theta[j] - e_phi[i] * e_phi[j]);
			h_cartesian[i][j] += h_cross * (e_theta[i] * e_phi[j] + e_phi[i] * e_theta[j]);
		}

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			cell.geom.tilde_gamma[i][j] += h_cartesian[i][j];

	double dh_plus_dt  = -omega * A / r * std::sin(phase) * envelope;
	double dh_cross_dt = -omega * A / r * std::cos(phase) * envelope;

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) {
			double dh_ij_dt = dh_plus_dt * (e_theta[i] * e_theta[j] - e_phi[i] * e_phi[j])
				+ dh_cross_dt * (e_theta[i] * e_phi[j] + e_phi[i] * e_theta[j]);
			cell.atilde.Atilde[i][j] += -0.5 / cell.gauge.alpha * dh_ij_dt;
		}
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

				double det_gamma = std::cbrt(
						cell.geom.gamma[0][0] * cell.geom.gamma[1][1] * cell.geom.gamma[2][2]
						- cell.geom.gamma[0][0] * cell.geom.gamma[1][2] * cell.geom.gamma[2][1]
						- cell.geom.gamma[1][1] * cell.geom.gamma[0][2] * cell.geom.gamma[2][0]
						- cell.geom.gamma[2][2] * cell.geom.gamma[0][1] * cell.geom.gamma[1][0]
						);
				cell.chi = 1.0 / det_gamma;

				for (int a = 0; a < 3; a++) {
					for (int b = 0; b < 3; b++) {
						cell.geom.tilde_gamma[a][b] = cell.chi * cell.geom.gamma[a][b];
					}
				}

				Matrix3x3 tilde_gamma_std;
				Matrix3x3 tilde_gamma_inv_std;

				for (int a = 0; a < 3; ++a)
					for (int b = 0; b < 3; ++b)
						tilde_gamma_std[a][b] = cell.geom.tilde_gamma[a][b];

				matrix.inverse_3x3(tilde_gamma_std, tilde_gamma_inv_std);

				for (int a = 0; a < 3; ++a)
					for (int b = 0; b < 3; ++b)
						cell.geom.tildgamma_inv[a][b] = tilde_gamma_inv_std[a][b];


				cell.gauge.alpha = 1.0 / std::sqrt(1.0 + 2.0 * H);
				cell.gauge.beta[0] = 2.0 * H * lx;
				cell.gauge.beta[1] = 2.0 * H * ly;
				cell.gauge.beta[2] = 2.0 * H * lz;

				for (int a_idx = 0; a_idx < 3; a_idx++) {
					for (int b_idx = 0; b_idx < 3; b_idx++) {
						cell.curv.K[a_idx][b_idx] = 0.0;
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
				/* 			printf("  %.5e", cell.atilde.Atilde[a][b]); */
				/* 		} */
				/* 		printf("\n"); */
				/* 	} */
				/* } */
				double trace_Atilde = 0.0;
				for (int a = 0; a < 3; ++a) {
					for (int b = 0; b < 3; ++b) {
						trace_Atilde += globalGrid[1][1][1].geom.tildgamma_inv[a][b] * globalGrid[1][1][1].atilde.Atilde[a][b];
					}
				}
                double Ktrace = 0.0;
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        Ktrace += globalGrid[i][j][k].geom.gamma_inv[a][b] * globalGrid[i][j][k].curv.K[a][b];
                    }
                }

				for (int a = 0; a < 3; a++) {
					for (int b = 0; b < 3; b++) {
						globalGrid[i][j][k].atilde.Atilde[a][b] = cell.chi * (globalGrid[i][j][k].curv.K[a][b] 
								- (1.0 / 3.0) * Ktrace * globalGrid[i][j][k].geom.gamma[a][b]);
					}
				}

            }
        }
	}
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
	printf("Gamma_inv_{ij} near (1,1,1) = \n");
	for (int p = 0; p < 3; p++) {
		printf("  ");
		for (int q = 0; q < 3; q++) {
			printf("%e ", globalGrid[1][1][1].geom.gamma_inv[p][q]);
		}
		printf("\n");
	}
	printf("Gamma_tile_{ij} near (1,1,1) = \n");
	for (int p = 0; p < 3; p++) {
		printf("  ");
		for (int q = 0; q < 3; q++) {
			printf("%e ", globalGrid[1][1][1].geom.tilde_gamma[p][q]);
		}
		printf("\n");
	}
	printf("Gamma_tilde_inv_{ij} near (1,1,1) = \n");
	for (int p = 0; p < 3; p++) {
		printf("  ");
		for (int q = 0; q < 3; q++) {
			printf("%e ", globalGrid[1][1][1].geom.tildgamma_inv[p][q]);
		}
		printf("\n");
	}

	printf("chi near (1,1,1) = %e\n", globalGrid[1][1][1].chi);

	printf("alpha_1_1_1 = %e\n", globalGrid[1][1][1].gauge.alpha);

	double test_radii[] = {0.5, 1.0, 2.0, 5.0, 10.0};
	for (double test_r : test_radii) {
		double H_test = M / test_r;  
		double alpha_test = 1.0 / std::sqrt(1.0 + 2.0 * H_test);
		printf("Eq plane r = %f : H = %e, alpha = %e (Schw approx)\n", 
				test_r, H_test, alpha_test);
	}
}
