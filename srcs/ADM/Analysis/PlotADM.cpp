
#include <Geodesics.h>

void export_gamma_slice(Grid &grid_obj, int j, double time) {
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

    double L = 9.0;
    double dx = (2.0 * L) / (NX - 1);
    double dz = (2.0 * L) / (NZ - 1);

    for (int i = 0; i < NX; ++i) {
        for (int k = 0; k < NZ; ++k) {
            double x = -L + i * dx;
            double z = -L + k * dz;

            const Grid::Cell2D &cell = grid_obj.getCell(i, j, k);

            file << x << "," << z << ","
                 << cell.geom.tilde_gamma[0][0] << "," << cell.geom.tilde_gamma[0][1] << "," << cell.geom.tilde_gamma[0][2] << ","
                 << cell.geom.tilde_gamma[1][0] << "," << cell.geom.tilde_gamma[1][1] << "," << cell.geom.tilde_gamma[1][2] << ","
                 << cell.geom.tilde_gamma[2][0] << "," << cell.geom.tilde_gamma[2][1] << "," << cell.geom.tilde_gamma[2][2] << "\n";
        }
    }

    file.close();
    std::cout << "[Export] Gamma slice saved to: " << filename.str() << std::endl;
}


void Grid::export_chi_slice(Grid &grid_obj, double time) {
    std::ofstream file;
    char filename[256];
    sprintf(filename, "Output/chi_slice_t%.3f.dat", time);
    file.open(filename);
	double L = 9.0;
    int z_mid = NZ / 2; // Coupe au plan z = 0

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double x = -L + i * DX;
            double y = -L + j * DY;
            double chi = globalGrid[i][j][z_mid].chi;
            file << x << " " << y << " " << chi << "\n";
        }
        file << "\n";
    }

    file.close();
}

void export_tilde_gamma_3D(Grid &grid_obj) {
    std::ofstream file("Output/tilde_gamma_full.vtk");
    file << "# vtk DataFile Version 2.0\n";
    file << "Conformal metric geom.tilde_gamma\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    file << "DIMENSIONS " << NX << " " << NY << " " << NZ << "\n";
    
    double x0 = -9.0;
    double y0 = -9.0;
    double z0 = -9.0;
    
    file << "ORIGIN " << x0 << " " << y0 << " " << z0 << "\n";

    double dx = 18.0 / (NX - 1);
    double dy = 18.0 / (NY - 1);
    double dz = 18.0 / (NZ - 1);
    
    file << "SPACING " << dx << " " << dy << " " << dz << "\n";

    file << "POINT_DATA " << (NX * NY * NZ) << "\n";
    file << "TENSORS geom.tilde_gamma float\n";

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
                file << cell.geom.tilde_gamma[0][0] << " " << cell.geom.tilde_gamma[0][1] << " " << cell.geom.tilde_gamma[0][2] << "\n";
                file << cell.geom.tilde_gamma[1][0] << " " << cell.geom.tilde_gamma[1][1] << " " << cell.geom.tilde_gamma[1][2] << "\n";
                file << cell.geom.tilde_gamma[2][0] << " " << cell.geom.tilde_gamma[2][1] << " " << cell.geom.tilde_gamma[2][2] << "\n\n";
            }
        }
    }

    // Ajout des composantes diagonales comme scalaires
    file << "SCALARS tilde_gamma_00 float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            for (int k = 0; k < NZ; k++)
                file << grid_obj.getCell(i, j, k).geom.tilde_gamma[0][0] << "\n";

    file << "SCALARS tilde_gamma_11 float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            for (int k = 0; k < NZ; k++)
                file << grid_obj.getCell(i, j, k).geom.tilde_gamma[1][1] << "\n";

    file << "SCALARS tilde_gamma_22 float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            for (int k = 0; k < NZ; k++)
                file << grid_obj.getCell(i, j, k).geom.tilde_gamma[2][2] << "\n";

    file.close();
    std::cout << "✅ Fichier VTK de la métrique conforme geom.tilde_gamma sauvegardé : tilde_gamma_full.vtk\n";
}

/* void export_constraints(std::string filename) { */
/*     std::ofstream file(filename); */
/*     file << "x,y,z,Hamiltonian,Momentum_x,Momentum_y,Momentum_z\n"; */
/*  */
/*     for (int i = 0; i < NX; i++) { */
/*         for (int j = 0; j < NY; j++) { */
/*             for (int k = 0; k < NZ; k++) { */
/*                 double x = -9.0 + i * (18.0 / (NX - 1)); */
/*                 double y = -9.0 + j * (18.0 / (NY - 1)); */
/*                 double z = -9.0 + k * (18.0 / (NZ - 1)); */
/*  */
/*                 double hamiltonian = hamiltonianGrid[i][j][k]; */
/*                 double momentum[3] = {0.0, 0.0, 0.0}; */
/*                 compute_constraints(*this, i, j, k, hamiltonian, momentum); */
/*  */
/*                 file << x << "," << y << "," << z << "," */
/*                      << hamiltonian << "," */
/*                      << momentum[0] << "," << momentum[1] << "," << momentum[2] << "\n"; */
/*             } */
/*         } */
/*     } */
/*  */
/*     file.close(); */
/*     std::cout << "✅ Contraintes exportées dans " << filename << std::endl; */
/* } */


void export_K_slice(Grid &grid_obj, int j) {
    std::ofstream file("Output/K_slice.csv");

    file << "x,z,K00,K01,K02,K10,K11,K12,K20,K21,K22\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            double x = -9.0 + i * (18.0 / (NX - 1));
            double z = -9.0 + k * (18.0 / (NZ - 1));

            Grid::Cell2D &cell = grid_obj.getCell(i, j, k);

            file << x << "," << z << ","
                 << cell.atilde.Atilde[0][0] << "," << cell.atilde.Atilde[0][1] << "," << cell.atilde.Atilde[0][2] << ","
                 << cell.atilde.Atilde[1][0] << "," << cell.atilde.Atilde[1][1] << "," << cell.atilde.Atilde[1][2] << ","
                 << cell.atilde.Atilde[2][0] << "," << cell.atilde.Atilde[2][1] << "," << cell.atilde.Atilde[2][2] << "\n";
        }
    }
    file.close();
    std::cout << "curv.K slice saved to K_slice.csv\n";
}

double r_e_plus(double theta, double a) {
    return M + sqrt(M * M - a * a * cos(theta) * cos(theta)); 
}

double r_e_minus(double theta, double a) {
    return M - sqrt(M * M - a * a * cos(theta) * cos(theta)); 
}

void export_K_3D(Grid &grid_obj) {
    std::ofstream file("Output/K_full.vtk");
    file << "# vtk DataFile Version 2.0\n";
    file << "curv.K extrinsic curvature\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    file << "DIMENSIONS " << NX << " " << NY << " " << NZ << "\n";

    double x0 = -9.0;
    double y0 = -9.0;
    double z0 = -9.0;
    double a = 0.9999;  
    file << "ORIGIN " << x0 << " " << y0 << " " << z0 << "\n";

    double dx = 18.0 / (NX - 1);
    double dy = 18.0 / (NY - 1);
    double dz = 18.0 / (NZ - 1);
    file << "SPACING " << dx << " " << dy << " " << dz << "\n";

    file << "POINT_DATA " << (NX * NY * NZ) << "\n";
    file << "TENSORS K float\n";

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
                file << cell.curv.K[0][0] << " " << cell.curv.K[0][1] << " " << cell.curv.K[0][2] << "\n";
                file << cell.curv.K[1][0] << " " << cell.curv.K[1][1] << " " << cell.curv.K[1][2] << "\n";
                file << cell.curv.K[2][0] << " " << cell.curv.K[2][1] << " " << cell.curv.K[2][2] << "\n\n";
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
    }    double r_H = M + sqrt(M * M - a * a); 

    file << "SCALARS Horizon float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = x0 + i * dx;
                double y = y0 + j * dy;
                double z = z0 + k * dz;
                double r = sqrt(x * x + y * y + z * z);
                file << (r < r_H ? 1.0 : 0.0) << "\n";
            }
        }
    }

    file << "SCALARS Ergosphere_plus float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = x0 + i * dx;
                double y = y0 + j * dy;
                double z = z0 + k * dz;
                double theta = atan2(sqrt(x * x + y * y), z);
                double r = sqrt(x * x + y * y + z * z);
                file << (r < r_e_plus(theta, a) ? 1.0 : 0.0) << "\n";
            }
        }
    }

    file << "SCALARS Ergosphere_minus float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = x0 + i * dx;
                double y = y0 + j * dy;
                double z = z0 + k * dz;
                double theta = atan2(sqrt(x * x + y * y), z);
                double r = sqrt(x * x + y * y + z * z);
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
            double x = -9.0 + i * (18.0 / (NX - 1));
            double z = -9.0 + k * (18.0 / (NZ - 1));

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
        for (int k = 0; k < NZ; k++) {
            double x = -9.0 + i * (18.0 / (NX - 1));
            double z = -9.0 + k * (18.0 / (NZ - 1));

            Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
            double d_alpha_dt, d_beta_dt[3];
            grid_obj.compute_gauge_derivatives(grid_obj, i, j, k, d_alpha_dt, d_beta_dt);

            file << x << "," << z << "," 
                 << cell.gauge.alpha << "," << cell.gauge.beta[0] << "," << cell.gauge.beta[1] << "," << cell.gauge.beta[2] << ","
                 << d_alpha_dt << "," << d_beta_dt[0] << "," << d_beta_dt[1] << "," << d_beta_dt[2] << "\n";
        }
    }

    file.close();
    std::cout << "Gauge slice saved to gauge_slice.csv\n";
}




void GridTensor::export_christoffel_slice(Grid &grid_obj, int j) {
    std::ofstream file("Output/christoffel_slice.csv");
	double L = 9.0;
    double x_min = -L, x_max = L;
    double y_min = -L, y_max = L;
    double z_min = -L, z_max = L;
    double dx = (x_max - x_min) / (NX - 1);
    double dy = (y_max - y_min) / (NY - 1);
    double dz = (z_max - z_min) / (NZ - 1);
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
            double x = x_min + i_idx * dx;
            double z = z_min + k_idx * dz;            
            double christof[3][3][3];
            compute_christoffel_3D(grid_obj, i_idx, j, k_idx, christof);
            
            file << x << "," << z;
            for (int i = 0; i < 3; i++) {
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        file << "," << christof[i][k][l];
                    }
                }
            }
            file << "\n";
        }
    }
    
    file.close();
    std::cout << "Christoffel slice saved to christoffel_slice.csv\n";
}




void Grid::export_fluid_slice(int j_slice) {
	double L = 2.0; 
    double x_min = -L, x_max = L;
    double y_min = -L, y_max = L;
    double z_min = -L, z_max = L;
    std::ofstream file("Output/fluid_slice.csv");
    if (!file.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier fluid_slice.csv" << std::endl;
        return;
    }

    file << "x,z,rho,p,matter.vx,vy,vz\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            double x = x_min + i * DX;
            double z = z_min + k * DZ;

            Cell2D &cell = globalGrid[i][j_slice][k];
            file << x << "," << z << ","
                 << cell.matter.rho << "," << cell.matter.p << ","
                 << cell.matter.vx << "," << cell.matter.vy << "," << cell.matter.vz
                 << "\n";
        }
    }

    file.close();
    std::cout << "✅ Export du fluide terminé dans fluid_slice.csv" << std::endl;
}



void Grid::export_energy_momentum_tensor_slice(int slice_y) {
    std::ofstream file("Output/T_energy_momentum.csv");
    if (!file.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier !\n";
        return;
    }

    file << "x,z";
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            file << ",T_" << a << b;
        }
    }
    file << "\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            Cell2D &cell = globalGrid[i][slice_y][k];

            file << i << "," << k;
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    file << "," << cell.matter.T[a][b];
                }
            }
            file << "\n";
        }
    }

    file.close();
    std::cout << "Exportation de T^{mu nu} terminée : T_energy_momentum.csv\n";
}

