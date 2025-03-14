#include <Geodesics.h>


double second_derivative_alpha(int dim_i, int dim_j,
                               int i, int j,
                               const std::vector<std::vector<double>> &alpha_grid,
                               double dr, double dtheta) {
    int Nr = alpha_grid.size();
    int Ntheta = (Nr > 0) ? alpha_grid[0].size() : 0;

    if (i <= 0 || i >= Nr-1 || j <= 0 || j >= Ntheta-1) {
        return 0.0;
    }
    
    if (dim_i == 0 && dim_j == 0) {
        return (alpha_grid[i+1][j] - 2.0 * alpha_grid[i][j] + alpha_grid[i-1][j]) / (dr * dr);
    }
    else if (dim_i == 1 && dim_j == 1) {
        return (alpha_grid[i][j+1] - 2.0 * alpha_grid[i][j] + alpha_grid[i][j-1]) / (dtheta * dtheta);
    }
    else if ((dim_i == 0 && dim_j == 1) || (dim_i == 1 && dim_j == 0)) {
        return (alpha_grid[i+1][j+1] - alpha_grid[i+1][j-1] - alpha_grid[i-1][j+1] + alpha_grid[i-1][j-1])
                / (4.0 * dr * dtheta);
		printf("alpha_grid[%d][%d] = %f\n", i, j, alpha_grid[i][j]);
    }
    return 0.0;
}

Matrix3x3 RK4_update_gamma(const Matrix3x3& gamma, const Matrix3x3& dgamma_dt, double dt) {
    Matrix3x3 result;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            result[i][j] = gamma[i][j] + dt * dgamma_dt[i][j];
    return result;
}

Matrix3x3 RK4_update_K(const Matrix3x3& K, const Matrix3x3& dK_dt, double dt) {
    Matrix3x3 result;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            result[i][j] = K[i][j] + dt * dK_dt[i][j];
    return result;
}

double compute_rho(double r) {
    double R = 5.0;
    double rho0 = 1e-3;
    return (r < R) ? rho0 : 0.0;
}



void matter_sources(double r, double& rho, Vector3& j, Matrix3x3& S_ij, double& S) {
    double R0 = 5.0;
    double sigma = 1.0;
    double rho0 = 1e-2; 
    rho = rho0 * exp(-pow((r - R0) / sigma, 2));

    j = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            S_ij[i][j] = 0.0;
        }
    }
    S = 0.0;
}


void evolveADM(Grid::Cell2D& cell, int i, int j, double dt, 
               const std::vector<std::vector<double>>& alpha_grid,
               double r_min, double theta_min, double dr, double dtheta,
               std::ofstream& file, int step) {

    Matrix3x3 R;
    Tensor3D Gamma3;
    Grid grid_obj;
    Matrix3x3 S_ij;
    double S;
    double rho;
    Vector3 j_i;

    double r = r_min + i * dr;
    double theta = theta_min + j * dtheta;

    matter_sources(r, rho, j_i, S_ij, S);
    cell.rho = rho;

    Vector3 X3D = { r, theta, 0.0 };
    grid_obj.calculate_christoffel_3D(X3D, Gamma3, cell.gamma, cell.gamma_inv);
    grid_obj.compute_ricci_3d(grid_obj, X3D, Gamma3, R);

    double K_trace = 0.0;
    for (int a = 0; a < 3; a++) {
        K_trace += cell.K[a][a];
    }

    Matrix3x3 dgamma_dt;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            dgamma_dt[a][b] = -2.0 * cell.alpha * cell.K[a][b];
        }
    }

    int i_index = i;
    int j_index = j;

    Matrix3x3 dK_dt;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double K2 = 0.0;
            for (int c = 0; c < 3; c++) {
                K2 += cell.K[a][c] * cell.K[c][b];
            }
            double D2alpha = 0.0;
            if (a < 2 && b < 2) {
                if (a == 0 && b == 0)
                    D2alpha = second_derivative_alpha(0, 0, i_index, j_index, alpha_grid, dr, dtheta);
                else if (a == 1 && b == 1)
                    D2alpha = second_derivative_alpha(1, 1, i_index, j_index, alpha_grid, dr, dtheta);
                else if ((a == 0 && b == 1) || (a == 1 && b == 0))
                    D2alpha = second_derivative_alpha(0, 1, i_index, j_index, alpha_grid, dr, dtheta);
            }
            double matter_term = 8.0 * M_PI * cell.alpha * (S_ij[a][b] - 0.5 * cell.gamma[a][b] * S);
            dK_dt[a][b] = -D2alpha + cell.alpha * (R[a][b] + K_trace * cell.K[a][b] - 2.0 * K2) + matter_term;
        }
    }

    cell.gamma = RK4_update_gamma(cell.gamma, dgamma_dt, dt);
    cell.K = RK4_update_K(cell.K, dK_dt, dt);

    file << step << "," << i << "," << j << "," << r << "," << theta;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            file << "," << dK_dt[a][b];
        }
    }
    file << "\n";
}
