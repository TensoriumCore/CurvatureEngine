#include <Geodesics.h>



void Grid::initializeFishboneMoncriefTorus(double r_in, double r_max, double rho_max, double l_torus, double Gamma) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                Cell2D &cell = globalGrid[i][j][k];

                double x = -9.0 + i * DX;
                double y = -9.0 + j * DY;
                double z = -9.0 + k * DZ;

                double r = sqrt(x * x + y * y);
                double R = sqrt(x * x + y * y + z * z);
                if (R < r_in || R > r_max) {
                    cell.matter.rho = 0.0;
                    cell.matter.p = 0.0;
                    cell.matter.vx = 0.0;
                    cell.matter.vy = 0.0;
                    cell.matter.vz = 0.0;
                    continue;
                }

                double theta = acos(z / (R + 1e-14));
                double sin_theta = sin(theta);
                double sin2_theta = sin_theta * sin_theta;

                double g_phiphi = r * r * sin2_theta;

                double u_phi_sq = (g_phiphi > 1e-14) ? (l_torus * l_torus / g_phiphi) : 0.0;

                double alpha = cell.gauge.alpha;
                double beta_sq = cell.gauge.beta[0] * cell.gauge.beta[0] +
                                 cell.gauge.beta[1] * cell.gauge.beta[1] +
                                 cell.gauge.beta[2] * cell.gauge.beta[2];

                double g_tt = -alpha * alpha + beta_sq;

                double W = -g_tt - u_phi_sq;

                if (W <= 0.0) {
                    cell.matter.rho = 0.0;
                    cell.matter.p = 0.0;
                    cell.matter.vx = 0.0;
                    cell.matter.vy = 0.0;
                    cell.matter.vz = 0.0;
                    continue;
                }

                double h = 1.0 + (Gamma / (Gamma - 1.0)) * (W - 1.0);
                double rho = rho_max * pow((h - 1.0) * (Gamma - 1.0) / Gamma, 1.0 / (Gamma - 1.0));
                double p = (Gamma - 1.0) * rho;

                cell.matter.rho = rho;
                cell.matter.p = p;

                double vphi = sqrt(u_phi_sq);
                cell.matter.vx = -y / r * vphi;
                cell.matter.vy =  x / r * vphi;
                cell.matter.vz = 0.0;
            }
        }
    }
}
