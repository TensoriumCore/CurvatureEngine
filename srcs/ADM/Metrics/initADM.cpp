#include "GridTensor.h"
#include <Geodesics.h>


void Grid::allocateGlobalGrid(){
    printf("Allocating global grid\n");

    globalGrid.resize(NX);
    for(int i = 0; i < NX; i++){
        globalGrid[i].resize(NY);
        for(int j = 0; j < NY; j++){
            globalGrid[i][j].resize(NZ);
        }
    }
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < NX; i++){
        for(int j = 0; j < NY; j++){
            for(int k = 0; k < NZ; k++){
                globalGrid[i][j][k] = Cell2D(); 
            }
        }
    }
}


void Grid::initializeData_Minkowski()
{
    printf("\n=== Initialisation Minkowski ===\n");

    double x_min = -128.0, x_max = 128.0;
    double y_min = -128.0, y_max = 128.0;
    double z_min = -128.0, z_max = 128.0;
    double dx = (x_max - x_min)/(NX-1);
    double dy = (y_max - y_min)/(NY-1);
    double dz = (z_max - z_min)/(NZ-1);

    for(int i=0; i<NX; i++)
    {
        double x = x_min + i*dx;
        for(int j=0; j<NY; j++)
        {
            double y = y_min + j*dy;
            for(int k=0; k<NZ; k++)
            {
                Cell2D &cell = globalGrid[i][j][k];

                cell.alpha = 1.0;
                cell.beta[0] = 0.0;
                cell.beta[1] = 0.0;
                cell.beta[2] = 0.0;

                for(int a=0; a<3; a++)
                {
                    for(int b=0; b<3; b++)
                    {
                        cell.gamma[a][b] = (a == b) ? 1.0 : 0.0;
                    }
                }

                for(int a=0; a<3; a++){
                    for(int b=0; b<3; b++){
                        cell.K[a][b] = 0.0;
                    }
                }
            }
        }
    }

    printf("Vérification Minkowski sur quelques points:\n");
    for(int test_i=0; test_i<3; test_i++)
    {
        for(int test_j=0; test_j<3; test_j++)
        {
            for(int test_k=0; test_k<3; test_k++)
            {
                Cell2D &cell = globalGrid[test_i][test_j][test_k];
                printf("Point (%d,%d,%d): alpha=%f, gamma=\n", test_i, test_j, test_k, cell.alpha);
                for(int a=0; a<3; a++)
                {
                    printf("  ");
                    for(int b=0; b<3; b++)
                    {
                        printf("%f ", cell.gamma[a][b]);
                    }
                    printf("\n");
                }
            }
        }
    }
    printf("=== Fin initialisation Minkowski ===\n");
}



void Grid::initializeKerrData(Grid &grid_obj) {
    double a = 0.998;   
    double L = 9.0;
    double x_min = -L, x_max = L;
    double y_min = -L, y_max = L;
    double z_min = -L, z_max = L;
    double dx = (x_max - x_min) / (NX - 1);
    double dy = (y_max - y_min) / (NY - 1);
    double dz = (z_max - z_min) / (NZ - 1);
    Matrix matrix;
	std::vector<std::array<double, 3>> horizonPoints;
    globalGrid.resize(NX);
    for (int i = 0; i < NX; i++) {
        globalGrid[i].resize(NY);
        for (int j = 0; j < NY; j++) {
            globalGrid[i][j].resize(NZ);
        }
    }

#pragma omp parallel for collapse(3) schedule(dynamic)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                double z = z_min + k * dz;
                double R2 = x * x + y * y + z * z;

				Cell2D &cell = globalGrid[i][j][k];

				double r2 = x*x + y*y + z*z;
				double r  = sqrt(r2);
				double r4 = r2*r2;
				double H  = (r < 1e-12) ? 0.0 
					: (M * r*r*r) / (r4 + a*a*z*z);

				double denom = (r2 + a*a);
				double lx = (r*x + a*y) / denom;
				double ly = (r*y - a*x) / denom;
				double lz = (r > 1e-12) ? (z / r) : 0.0;

				double lxlx = lx * lx;
				double lyly = ly * ly;
				double lzlz = lz * lz;

				cell.gamma[0][0] = 1.0 + 2.0*H*lxlx; 
				cell.gamma[0][1] = 2.0*H*lx*ly;      
				cell.gamma[0][2] = 2.0*H*lx*lz;
				cell.gamma[1][0] = cell.gamma[0][1]; 
				cell.gamma[1][1] = 1.0 + 2.0*H*lyly; 
				cell.gamma[1][2] = 2.0*H*ly*lz;
				cell.gamma[2][0] = cell.gamma[0][2];
				cell.gamma[2][1] = cell.gamma[1][2];
				cell.gamma[2][2] = 1.0 + 2.0*H*lzlz;


                cell.alpha = 1.0 / sqrt(1.0 + 2 * H);

                matrix.inverse_3x3(cell.gamma, cell.gamma_inv);
				cell.beta[0] = 2.0 * H * lx;
				cell.beta[1] = 2.0 * H * ly;
				cell.beta[2] = 2.0 * H * lz;

				cell.alpha = 1.0 / sqrt(1.0 + 2.0 * H);

				for (int a_idx = 0; a_idx < 3; a_idx++) {
					for (int b_idx = 0; b_idx < 3; b_idx++) {
						cell.K[a_idx][b_idx] = 0.0;
					}
				}
				#ifndef PFLUID
				
					double r_cart = sqrt(x * x + y * y + z * z);
					cell.rho = exp(-r_cart * r_cart / 2.0);
					cell.p = 0.3 * cell.rho + 0.5 * cell.rho * cell.rho;

					double vr = 2.0;  
					if (r_cart > 1e-6) {
						cell.vx = -vr * y / r_cart;
						cell.vy = vr * x / r_cart;
						cell.vz = vr * z / r_cart * (1.0 - 1.0 / r_cart);
					} else {
						cell.vx = cell.vy = 0.0;
					}
					cell.vz = 0.0;

					double r_horizon = M + sqrt(M * M - a * a);
					double epsilon = 1e-2;

					double diff = fabs(r - r_horizon);
					if (diff < epsilon) {
						horizonPoints.push_back({x, y, z});
					}
					/* if (r < r_horizon) { */
					/* 	cell.vx = cell.vy = cell.vz = 0.0; */
					/* 	cell.rho = 0.0; */
					/* 	cell.p = 0.0; */
					/* } */

				#endif

                cell.beta[0] = 2 * H * lx;
                cell.beta[1] = 2 * H * ly;
                cell.beta[2] = 2 * H * lz;
				double norm_beta = sqrt(cell.beta[0] * cell.beta[0] + cell.beta[1] * cell.beta[1] + cell.beta[2] * cell.beta[2]);
				/* printf("norm_beta = %e\n", norm_beta); */
				if (norm_beta > 1.0) {
					double scale = 1.0 / norm_beta;
					cell.beta[0] *= scale;
					cell.beta[1] *= scale;
					cell.beta[2] *= scale;

				}
				/* if (i == NX/2 && j == NY/2 && k == NZ/2) { */
				/* 	printf("beta[0](+dx) = %e, beta[0](-dx) = %e\n", */
				/* 			grid_obj.getCell(i+1, j, k).beta[0], grid_obj.getCell(i-1, j, k).beta[0]); */
				/* 	printf("beta[1](+dy) = %e, beta[1](-dy) = %e\n", */
				/* 			grid_obj.getCell(i, j+1, k).beta[1], grid_obj.getCell(i, j-1, k).beta[1]); */
				/* 	printf("beta[2](+dz) = %e, beta[2](-dz) = %e\n", */
				/* 			grid_obj.getCell(i, j, k+1).beta[2], grid_obj.getCell(i, j, k-1).beta[2]); */
				/* } */
				if (j == NY / 2 && k == NZ / 2) {
                    printf("gamma[0][0] à (i=%d, j=%d, k=%d) = %e\n", i, j, k, cell.gamma[0][0]);
                }
				if (i == NX / 2 && j == NY / 2 && k == NZ / 2) {
					printf("H = %e at r = %e (x=%e, y=%e, z=%e)\n", H, r, x, y, z);
				}
			}
		}
	}
	GridTensor gridtensor;
	for (int i = 1; i < NX - 1; i++) {
		for (int j = 1; j < NY - 1; j++) {
			for (int k = 1; k < NZ - 1; k++) {
				gridtensor.compute_extrinsic_curvature(grid_obj, i, j, k, dx, dy, dz);
			}
		}
	}
	printf("K = \n");
	for (int i = 0; i < 3; i++) {
		printf("  ");
		for (int j = 0; j < 3; j++) {
			printf("%e ", globalGrid[1][1][1].K[i][j]);
		}
		printf("\n");
	}
	
	printf("alpha_1_1_1 = %e\n", globalGrid[1][1][1].alpha);
    double test_radii[] = {0.5, 1.0, 2.0, 5.0, 10.0};
    for (double test_r : test_radii) {
        double H_test = M / test_r;
        double alpha_test = 1.0 / sqrt(1.0 + 2 * H_test);
        printf("Eq plane r = %f : H = %e, alpha = %e\n", test_r, H_test, alpha_test);
    }
	{
        std::ofstream ofs("horizon.csv");
        ofs << "x,y,z\n";
        for (auto &pt : horizonPoints) {
            ofs << pt[0] << "," << pt[1] << "," << pt[2] << "\n";
        }
        ofs.close();
    }
}


void Grid::initializeData() {
    double L = 4.0; 
    double x_min = -L, x_max = L;
    double y_min = -L, y_max = L;
    double z_min = -L, z_max = L;
    double dx = (x_max - x_min) / (NX - 1);
    double dy = (y_max - y_min) / (NY - 1);
    double dz = (z_max - z_min) / (NZ - 1);
    Matrix matrix;

    globalGrid.resize(NX);
    for (int i = 0; i < NX; i++) {
        globalGrid[i].resize(NY);
        for (int j = 0; j < NY; j++) {
            globalGrid[i][j].resize(NZ);
        }
    }

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                double z = z_min + k * dz;
                double r = sqrt(x * x + y * y + z * z);

                Cell2D cell;
                double Phi = 1.0 + 0.5 * M / r;
                cell.alpha = (1.0 - M / (2 * r)) / (1.0 + M / (2 * r));
                
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        cell.gamma[a][b] = (a == b) ? pow(Phi, 4) : 0.0;
                    }
                }

                matrix.inverse_3x3(cell.gamma, cell.gamma_inv);
                
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        cell.K[a][b] = 0.0;
                    }
                }

                cell.rho = exp(-r * r / 2.0);
                cell.p = 0.3 * cell.rho + 0.5 * cell.rho * cell.rho;

                double vr = 2.2;
                cell.vx = vr * x / r;
                cell.vy = vr * y / r;
                cell.vz = vr * z / r;

                if (r < 2 * M) {
                    cell.vx = 0.0;
                    cell.vy = 0.0;
                    cell.vz = 0.0;
                    cell.rho = 0.0;
                    cell.p = 0.0;
                }

                globalGrid[i][j][k] = cell;

                if (j == NY / 2 && k == NZ / 2) { 
                    printf("gamma[0][0] at (i=%d, j=%d, k=%d) = %e\n", i, j, k, globalGrid[i][j][k].gamma[0][0]);
                }
            }
        }
    }
	double test_radii[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 128.0};
	for (double test_r : test_radii) {
		double Phi_test = 1.0 + 0.5 * M / test_r;
		double alpha_test = (1.0 - M / (2 * test_r)) / (1.0 + M / (2 * test_r));
		printf("r = %f : Phi^4 = %e, alpha = %e\n", test_r, pow(Phi_test, 4), alpha_test);
	}

	for (int iCell = 0; iCell < 3; iCell++) {
		for (int jCell = 0; jCell < 3; jCell++) {
			for (int kCell = 0; kCell < 3; kCell++) {
				Cell2D &cell = globalGrid[iCell][jCell][kCell];
				for (int a = 0; a < 3; a++) {
					for (int b = 0; b < 3; b++) {
						if (a != b && fabs(cell.gamma[a][b]) > 1e-12) {
							printf("⚠️ γ[%d][%d] non diagonale en (%d,%d,%d) : %e\n", 
									a, b, iCell, jCell, kCell, cell.gamma[a][b]);
						}
					}
				}
			}
		}
	}
	int i_far = NX - 1;  
	int j_center = NY / 2;
	int k_center = NZ / 2;
	Cell2D &cell_far = globalGrid[i_far][j_center][k_center];

	std::cout << "Test 1 (r -> ∞):\n";
	std::cout << "Alpha = " << cell_far.alpha << "\n";
	std::cout << "Gamma_xx = " << cell_far.gamma[0][0] << "\n";

	double rho_horizon = M / 2.0;  
	int i_horizon = static_cast<int>((rho_horizon - x_min) / dx);
	Cell2D &cell_horizon = globalGrid[i_horizon][j_center][k_center];

	std::cout << "\nTest 2 (r = M/2):\n";
	std::cout << "Alpha = " << cell_horizon.alpha << "\n";
	std::cout << "Position x = " << (x_min + i_horizon*dx) << "\n";
}




