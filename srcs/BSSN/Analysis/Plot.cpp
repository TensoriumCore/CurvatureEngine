
#include <Geodesics.h>

void export_gamma_slice(Grid &grid_obj, int j, float time) {
    std::ostringstream filename;
    filename << "Output/gamma_slice_t" << std::fixed << std::setprecision(3) << time << ".csv";

    std::ofstream file(filename.str());
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename.str() << std::endl;
        return;
    }

    file << "x,z,"
         << "gamma_00,gamma_01,gamma_02,"
         << "gamma_10,gamma_11,gamma_12,"
         << "gamma_20,gamma_21,gamma_22\n";

    float L = 9.0;
    float dx = (2.0 * L) / (NX - 1);
    float dz = (2.0 * L) / (NZ - 1);

    for (int i = 0; i < NX; ++i) {
        for (int k = 0; k < NZ; ++k) {
            float x = -L + i * dx;
            float z = -L + k * dz;

            const Grid::Cell2D &cell = grid_obj.getCell(i, j, k);

            file << x << "," << z << ","
                 << cell.geom.dt_tilde_gamma[0][0] << "," << cell.geom.dt_tilde_gamma[0][1] << "," << cell.geom.dt_tilde_gamma[0][2] << ","
                 << cell.geom.dt_tilde_gamma[1][0] << "," << cell.geom.dt_tilde_gamma[1][1] << "," << cell.geom.dt_tilde_gamma[1][2] << ","
                 << cell.geom.dt_tilde_gamma[2][0] << "," << cell.geom.dt_tilde_gamma[2][1] << "," << cell.geom.dt_tilde_gamma[2][2] << "\n";
        }
    }

    file.close();
    std::cout << "[Export] Gamma slice saved to: " << filename.str() << std::endl;
}

void Grid::appendConstraintL2ToCSV(const std::string& filename, float time) const {
    float sum_H = 0.0;
    float sum_Mx = 0.0, sum_My = 0.0, sum_Mz = 0.0;
    int N = 0;

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int k = 1; k < NZ - 1; ++k) {
                const Cell2D &cell = globalGrid[i][j][k];
                sum_H  += cell.matter.hamiltonian * cell.matter.hamiltonian;
                sum_Mx += cell.matter.momentum[0] * cell.matter.momentum[0];
                sum_My += cell.matter.momentum[1] * cell.matter.momentum[1];
                sum_Mz += cell.matter.momentum[2] * cell.matter.momentum[2];
                ++N;
            }
        }
    }

    float L2_H  = std::sqrt(sum_H / N);
    float L2_Mx = std::sqrt(sum_Mx / N);
    float L2_My = std::sqrt(sum_My / N);
    float L2_Mz = std::sqrt(sum_Mz / N);

    std::ofstream file;
    bool exists = std::ifstream(filename).good();
    file.open(filename.c_str(), std::ios::app);

    if (!exists) {
        file << "time,hamiltonian_L2,momentum_x_L2,momentum_y_L2,momentum_z_L2\n";
    }

    file << time << ","
         << L2_H << ","
         << L2_Mx << ","
         << L2_My << ","
         << L2_Mz << "\n";

    file.close();
}



void Grid::export_Atildedt_slide(Grid &grid_obj, float time) {
	std::ofstream file;
	char filename[256];
	sprintf(filename, "Output/Atildedt_slice_t%.3f.dat", time);
	file.open(filename);
	float L = 19.0;
	int z_mid = NY / 2;

	for (int i = 0; i < NX; ++i) {
		for (int j = 0; j < NY; ++j) {
			float x = -L + i * DX;
			float y = -L + j * DY;
			float Atildedt = globalGrid[i][j][z_mid].atilde.dt_Atilde[0][0];
			file << x << " " << y << " " << Atildedt << "\n";
		}
		file << "\n";
	}
}


void export_K_slice(Grid &grid_obj, int j) {
    std::ofstream file("Output/K_slice.csv");

    file << "x,z,K00,K01,K02,K10,K11,K12,K20,K21,K22\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            float x = -12.0 + i * (9.0 / (NX - 1));
            float z = -12.0 + k * (9.0 / (NZ - 1));

            Grid::Cell2D &cell = grid_obj.getCell(i, j, k);

            file << x << "," << z << ","
                 << cell.atilde.dt_Atilde[0][0] << "," << cell.atilde.dt_Atilde[0][1] << "," << cell.atilde.dt_Atilde[0][2] << ","
                 << cell.atilde.dt_Atilde[1][0] << "," << cell.atilde.dt_Atilde[1][1] << "," << cell.atilde.dt_Atilde[1][2] << ","
                 << cell.atilde.dt_Atilde[2][0] << "," << cell.atilde.dt_Atilde[2][1] << "," << cell.atilde.dt_Atilde[2][2] << "\n";
        }
    }
    file.close();
    std::cout << "curv.K slice saved to K_slice.csv\n";
}

float r_e_plus(float theta, float a) {
    return M + sqrt(M * M - a * a * cos(theta) * cos(theta)); 
}

float r_e_minus(float theta, float a) {
    return M - sqrt(M * M - a * a * cos(theta) * cos(theta)); 
}

void export_K_3D(Grid &grid_obj) {
    std::ofstream file("Output/K_full.vtk");
    file << "# vtk DataFile Version 2.0\n";
    file << "curv.K extrinsic curvature\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    file << "DIMENSIONS " << NX << " " << NY << " " << NZ << "\n";

    float x0 = -9.0;
    float y0 = -9.0;
    float z0 = -9.0;
    float a = 0.9999;  
    file << "ORIGIN " << x0 << " " << y0 << " " << z0 << "\n";

    float dx = 18.0 / (NX - 1);
    float dy = 18.0 / (NY - 1);
    float dz = 18.0 / (NZ - 1);
    file << "SPACING " << dx << " " << dy << " " << dz << "\n";

    file << "POINT_DATA " << (NX * NY * NZ) << "\n";
    file << "TENSORS K float\n";

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
                file << cell.atilde.Atilde[0][0] << " " << cell.atilde.Atilde[0][1] << " " << cell.atilde.Atilde[0][2] << "\n";
                file << cell.atilde.Atilde[1][0] << " " << cell.atilde.Atilde[1][1] << " " << cell.atilde.Atilde[1][2] << "\n";
                file << cell.atilde.Atilde[2][0] << " " << cell.atilde.Atilde[2][1] << " " << cell.atilde.Atilde[2][2] << "\n\n";
            }
        }
    }
    file << "TENSORS dKt float\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
                file << cell.atilde.dt_Atilde[0][0] << " " << cell.atilde.dt_Atilde[0][1] << " " << cell.atilde.dt_Atilde[0][2] << "\n";
                file << cell.atilde.dt_Atilde[1][0] << " " << cell.atilde.dt_Atilde[1][1] << " " << cell.atilde.dt_Atilde[1][2] << "\n";
                file << cell.atilde.dt_Atilde[2][0] << " " << cell.atilde.dt_Atilde[2][1] << " " << cell.atilde.dt_Atilde[2][2] << "\n\n";
            }
        }
    }    float r_H = M + sqrt(M * M - a * a); 

    file << "SCALARS Horizon float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                float x = x0 + i * dx;
                float y = y0 + j * dy;
                float z = z0 + k * dz;
                float r = sqrt(x * x + y * y + z * z);
                file << (r < r_H ? 1.0 : 0.0) << "\n";
            }
        }
    }

    file << "SCALARS Ergosphere_plus float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                float x = x0 + i * dx;
                float y = y0 + j * dy;
                float z = z0 + k * dz;
                float theta = atan2(sqrt(x * x + y * y), z);
                float r = sqrt(x * x + y * y + z * z);
                file << (r < r_e_plus(theta, a) ? 1.0 : 0.0) << "\n";
            }
        }
    }

    file << "SCALARS Ergosphere_minus float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                float x = x0 + i * dx;
                float y = y0 + j * dy;
                float z = z0 + k * dz;
                float theta = atan2(sqrt(x * x + y * y), z);
                float r = sqrt(x * x + y * y + z * z);
                file << (r < r_e_minus(theta, a) ? 1.0 : 0.0) << "\n";
            }
        }
    }
	
	file << "SCALARS fluid float 1\n";
	file << "LOOKUP_TABLE default\n";
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NY; j++) {
			for (int k = 0; k < NZ; k++) {
				Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
				file << (cell.matter.rho > 1e-10 ? 1.0 : 0.0) << "\n";
			}
		}
	}
	
	file << "SCALARS fluid_velocity float 1\n";
	file << "LOOKUP_TABLE default\n";
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NY; j++) {
			for (int k = 0; k < NZ; k++) {
				Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
				file << sqrt(cell.matter.vx * cell.matter.vx + cell.matter.vy * \
						cell.matter.vy + cell.matter.vz * cell.matter.vz) << "\n";
			}
		}
	}

	file << "SCALARS alpha float 1\n";
	file << "LOOKUP_TABLE default\n";
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NY; j++) {
			for (int k = 0; k < NZ; k++) {
				auto &cell = grid_obj.getCell(i, j, k);
				file << cell.gauge.alpha << "\n";
			}
		}
	}
	file << "VECTORS shift float\n";
for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
        for (int k = 0; k < NZ; k++) {
            auto &cell = grid_obj.getCell(i, j, k);
            file << cell.gauge.beta[0] << " " << cell.gauge.beta[1] << " " << cell.gauge.beta[2] << "\n";
        }
    }
}

    file.close();
    std::cout << "curv.K 3D VTK file with Kerr surfaces saved to K_full.vtk\n";
}


void export_alpha_slice(Grid &grid_obj, int j) {
    std::ofstream file("Output/alpha_slice.csv");

    file << "x,z,alpha\n";
    for(int i = 0; i < NX; i++) {
        for(int k = 0; k < NZ; k++) {
            float x = -9.0 + i * (18.0 / (NX - 1));
            float z = -9.0 + k * (18.0 / (NZ - 1));

            Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
            file << x << "," << z << "," << cell.gauge.alpha << "\n";
        }
    }
    file.close();
    std::cout << "Alpha slice saved to alpha_slice.csv\n";
}



void export_gauge_slice(Grid &grid_obj, int j) {
    std::ofstream file("Output/gauge_slice.csv");
    file << "x,z,alpha,beta0,beta1,beta2,d_alpha_dt,d_beta0_dt,d_beta1_dt,d_beta2_dt\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NY; k++) {
            float x = -9.0 + i * (18.0 / (NX - 1));
            float z = -9.0 + k * (18.0 / (NZ - 1));

            Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
            float d_alpha_dt, d_beta_dt[3];
            grid_obj.compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);

            file << x << "," << z << "," 
                 << cell.gauge.alpha << "," << cell.gauge.beta[0] << "," << cell.gauge.beta[1] << "," << cell.gauge.beta[2] << ","
                 << cell.gauge.dt_alpha << "," << cell.gauge.dt_beta[0] << "," << cell.gauge.dt_beta[1] << "," << cell.gauge.dt_beta[2] << "\n";
        }
    }

    file.close();
    std::cout << "Gauge slice saved to gauge_slice.csv\n";
}


void export_gauge_slice_anim(Grid &grid_obj, int j, float time) {
    char name[256];
    sprintf(name, "Output/gauge_slice_t%.3f.csv", time);
    std::ofstream file(name);
    file << "x,z,alpha,beta0,beta1,beta2,d_alpha_dt,d_beta0_dt,d_beta1_dt,d_beta2_dt\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            float x = -9.0f + i * (18.0f / (NX - 1));
            float z = -9.0f + k * (18.0f / (NZ - 1));

            Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
            float d_alpha_dt, d_beta_dt[3];
            grid_obj.compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);

            file << x << "," << z << ","
                 << cell.gauge.alpha << ","
                 << cell.gauge.beta[0] << "," << cell.gauge.beta[1] << "," << cell.gauge.beta[2] << ","
                 << d_alpha_dt << "," << cell.gauge.beta[0] << "," << cell.gauge.beta[1] << "," << cell.gauge.beta[2] << "\n";
        }
    }
    file.close();
    std::cout << "[Export] Gauge slice saved to: " << name << "\n";
}

void export_gauge_slice_xy(Grid &grid, int k, float time) {
    char name[256];
    sprintf(name, "Output/gauge_xy_t%.3f.csv", time);
    std::ofstream f(name);
    f << "x,y,alpha,beta0,beta1\n";

    float L=9.f, dx=2*L/(NX-1), dy=2*L/(NY-1);
    for (int i=0;i<NX;++i)
      for (int j=0;j<NY;++j) {
        float x=-L+i*dx, y=-L+j*dy;
        const auto &c = grid.getCell(i,j,k);
        f << x << "," << y << "," << c.gauge.alpha << ","
          << c.gauge.dt_beta[0] << "," << c.gauge.dt_beta[1] << "\n";
      }
    printf("[Export] Gauge XY slice saved: %s\n", name);
}

void GridTensor::export_christoffel_slice(Grid &grid_obj, int j) {
    std::ofstream file("Output/christoffel_slice.csv");
	float L = 9.0;
    float x_min = -L, x_max = L;
    float y_min = -L, y_max = L;
    float z_min = -L, z_max = L;
    float dx = (x_max - x_min) / (NX - 1);
    float dy = (y_max - y_min) / (NY - 1);
    float dz = (z_max - z_min) / (NZ - 1);
    file << "x,z";
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                file << ",Gamma" << i << k << l;
            }
        }
    }
    file << "\n";
    
    for (int i_idx = 1; i_idx < NX-1; i_idx++) {
        for (int k_idx = 1; k_idx < NZ-1; k_idx++) {
            float x = x_min + i_idx * dx;
            float z = z_min + k_idx * dz;            
            float christof[3][3][3];
			Grid::Cell2D &cell = grid_obj.getCell(i_idx, j, k_idx);

            file << x << "," << z;
            for (int i = 0; i < 3; i++) {
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        file << "," << cell.conn.Christoffel[i][k][l];
                    }
                }
            }
            file << "\n";
        }
    }
    
    file.close();
    std::cout << "Christoffel slice saved to christoffel_slice.csv\n";
}


