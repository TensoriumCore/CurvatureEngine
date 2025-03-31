#include <Geodesics.h>
#include <tuple>


static inline double kerrSchildRadius(double x, double y, double z, double a)
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

void Grid::initializeBinaryKerrData(Grid &grid_obj) {
    double m1 = 1.0, a1 = 0.935;
    double m2 = 1.0, a2 = 0.935;

    double x1 = 0.0, y1 = -3.0, z1 = 0.0; 
    double x2 = 0.0, y2 = 3.0, z2 = 0.0; 

    double L = 24.0;
    double x_min = -L, x_max = L;
    double y_min = -L, y_max = L;
    double z_min = -L, z_max = L;

    double dxBH = x2 - x1;
    double dyBH = y2 - y1;
    double dzBH = z2 - z1;
    double r12  = std::sqrt(dxBH*dxBH + dyBH*dyBH + dzBH*dzBH);
    double eta  = m1 * m2 / (m1 + m2);
    double v_orb = std::sqrt((m1 + m2)/r12) * (1.0 + (3.0 - eta)/r12);

    double dx = (x_max - x_min) / (NX - 1);
    double dy = (y_max - y_min) / (NY - 1);
    double dz = (z_max - z_min) / (NZ - 1);

    GridTensor gridtensor;
    Matrix matrix;
    globalGrid.resize(NX, std::vector<std::vector<Cell2D>>(NY, std::vector<Cell2D>(NZ)));

    
	auto lorentz_boost_lmu = [](double beta, double lt, double lx, double ly, double lz) {
		const double gamma = 1.0 / sqrt(1.0 - beta * beta);
		double lt_prime = gamma * (lt + beta * ly);
		double ly_prime = gamma * (ly + beta * lt);
		return std::make_tuple(lt_prime, lx, ly_prime, lz);
	};

    #pragma omp parallel for collapse(3) schedule(dynamic)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                double z = z_min + k * dz;

                double dx1 = x - x1, dy1 = y - y1, dz1 = z - z1;
                double rKS1 = kerrSchildRadius(dx1, dy1, dz1, a1);
                double cosTheta1 = (rKS1 > 1e-14) ? dz1 / rKS1 : 0.0;
                double denomKS1 = rKS1*rKS1 + a1*a1*cosTheta1*cosTheta1;
                double H1 = (rKS1 > 1e-14 && denomKS1 > 1e-14) ? (m1 * rKS1) / denomKS1 : 0.0;
                double denomVec1 = rKS1*rKS1 + a1*a1;
                double lx1 = (denomVec1 > 1e-14) ? (rKS1 * dx1 + a1 * dy1) / denomVec1 : 0.0;
                double ly1 = (denomVec1 > 1e-14) ? (rKS1 * dy1 - a1 * dx1) / denomVec1 : 0.0;
                double lz1 = (rKS1 > 1e-14) ? dz1 / rKS1 : 0.0;

				auto [lt1, lx1b, ly1b, lz1b] = lorentz_boost_lmu(+v_orb, 1.0, lx1, ly1, lz1);

                double dx2 = x - x2, dy2 = y - y2, dz2 = z - z2;
                double rKS2 = kerrSchildRadius(dx2, dy2, dz2, a2);
                double cosTheta2 = (rKS2 > 1e-14) ? dz2 / rKS2 : 0.0;
                double denomKS2 = rKS2*rKS2 + a2*a2*cosTheta2*cosTheta2;
                double H2 = (rKS2 > 1e-14 && denomKS2 > 1e-14) ? (m2 * rKS2) / denomKS2 : 0.0;
                double denomVec2 = rKS2*rKS2 + a2*a2;
                double lx2 = (denomVec2 > 1e-14) ? (rKS2 * dx2 + a2 * dy2) / denomVec2 : 0.0;
                double ly2 = (denomVec2 > 1e-14) ? (rKS2 * dy2 - a2 * dx2) / denomVec2 : 0.0;
                double lz2 = (rKS2 > 1e-14) ? dz2 / rKS2 : 0.0;
				auto [lt2b, lx2b, ly2b, lz2b] = lorentz_boost_lmu(-v_orb, 1.0, lx2, ly2, lz2);

                double H = H1 + H2;
                double lx = lx1b + lx2b;
                double ly = ly1b + ly2b;
                double lz = lz1b + lz2b;
                double norm_l = std::sqrt(lx*lx + ly*ly + lz*lz);
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
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        cell.geom.tilde_gamma[a][b] = cell.geom.gamma[a][b];

                Matrix3x3 tg_std, inv_tg_std;
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        tg_std[a][b] = cell.geom.tilde_gamma[a][b];
                matrix.inverse_3x3(tg_std, inv_tg_std);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        cell.geom.tildgamma_inv[a][b] = inv_tg_std[a][b];

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


	BSSNevolve bssn;
	for (int i = 1; i < NX - 1; i++) {
		for (int j = 1; j < NY - 1; j++) {
			for (int k = 1; k < NZ - 1; k++) {
				for (int a = 0; a < 3; ++a) {
					for (int b = 0; b < 3; ++b) {
						globalGrid[i][j][k].dgt[a][b] = 0.0;  // dt_tilde_gamma = 0
					}
				}		
				gridtensor.compute_extrinsic_curvature(*this, i, j, k, dx, dy, dz);
				gridtensor.compute_Atilde(*this, i, j, k);

				Cell2D &cell = globalGrid[i][j][k];
				double Ktrace = 0.0;
				for (int a = 0; a < 3; a++)
					for (int b = 0; b < 3; b++)
						Ktrace += cell.geom.tildgamma_inv[a][b] * cell.curv.K[a][b];
				for (int a = 0; a < 3; a++)
					for (int b = 0; b < 3; b++)
						cell.atilde.Atilde[a][b] = cell.chi * (cell.curv.K[a][b] - (1.0 / 3.0) * Ktrace * cell.geom.tilde_gamma[a][b]);
			}
		}
	}


	printf("Finished computing Atilde and K\n");

	int max_iter = 5000;
	double tol = 1e-8;
	solve_lichnerowicz(max_iter, tol, dx, dy, dz);
	printf("Chi = %f\n", globalGrid[1][1][1].chi);
	for (int i = 1; i < NX - 1; i++) {
		for (int j = 1; j < NY - 1; j++) {
			for (int k = 1; k < NZ - 1; k++) {
				Cell2D &cell = globalGrid[i][j][k];
				for (int a = 0; a < 3; a++)
					for (int b = 0; b < 3; b++)
						cell.geom.tilde_gamma[a][b] = cell.chi * cell.geom.gamma[a][b];
				Matrix3x3 tg_std, inv_tg_std;
				for (int a = 0; a < 3; ++a)
					for (int b = 0; b < 3; ++b)
						tg_std[a][b] = cell.geom.tilde_gamma[a][b];
				matrix.inverse_3x3(tg_std, inv_tg_std);
				for (int a = 0; a < 3; ++a)
					for (int b = 0; b < 3; ++b)
						cell.geom.tildgamma_inv[a][b] = inv_tg_std[a][b];
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
    printf("chi near (1,1,1) = %e\n", globalGrid[1][1][1].chi);
    printf("alpha_1_1_1 = %e\n", globalGrid[1][1][1].gauge.alpha);
    
    double test_radii[] = {0.5, 1.0, 2.0, 5.0, 10.0};
    for (double test_r : test_radii) {
        double H_test = M / test_r;  
        double alpha_test = 1.0 / std::sqrt(1.0 + 2.0 * H_test);
        printf("Eq plane r = %f : H = %e, alpha = %e (Schw approx)\n", test_r, H_test, alpha_test);
    }
}
