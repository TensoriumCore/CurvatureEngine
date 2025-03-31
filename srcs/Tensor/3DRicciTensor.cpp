#include <Geodesics.h>



void Grid::calculate_ricci_3d_from_riemann(const Riemann3D& Riemann, Matrix3x3& Ricci) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += Riemann[k][i][k][j];
            }
            Ricci[i][j] = sum;
        }
    }
    printf("\nRicci tensor:\n");
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            printf("%12.6f\t", Ricci[mu][nu]);
        }
        printf("\n");
    }
    

}




void compute_partial_christoffel_3D(
    const Vector3& X,   
    int m,             
    Tensor3D& dGamma,  
    double delta
) {
    Vector3 Xp = X;
    Vector3 Xm = X;
    Grid grid_obj; 
    Xp[m] += delta;
    Xm[m] -= delta;
    
    Tensor3D Gamma_p{}; 
    Tensor3D Gamma_m{};
    Matrix3x3 gamma{};   
    Matrix3x3 gamma_inv{};

    grid_obj.calculate_christoffel_3D(Xp, Gamma_p, gamma, gamma_inv);
    grid_obj.calculate_christoffel_3D(Xm, Gamma_m, gamma, gamma_inv);

    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dGamma[k][i][j] = (Gamma_p[k][i][j] - Gamma_m[k][i][j]) / (2.0 * delta);
            }
        }
    }
}


void Grid::compute_ricci_3d(
    const Vector3& X,       
    const Tensor3D& Gamma3, 
    Matrix3x3& R3   
) {
    for (auto &row : R3) {
        row.fill(0.0);
    }

    static Tensor4D partialGamma{}; 

    double delta = 1e-5;

    for (int m = 0; m < 3; m++) {
        std::array<std::array<std::array<double, 3>, 3>, 3> dG{};
        compute_partial_christoffel_3D(X, m, dG, delta);
        
        for (int k = 0; k < 3; k++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    partialGamma[m][k][i][j] = dG[k][i][j];
                }
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;
            for (int k = 0; k < 3; k++) {
                term1 += partialGamma[k][k][i][j];
                term2 += partialGamma[j][k][i][k];
            }
            for (int k = 0; k < 3; k++) {
                for (int m = 0; m < 3; m++) {
                    term3 += Gamma3[k][i][j] * Gamma3[m][k][m]; 
                }
            }
            for (int m = 0; m < 3; m++) {
                for (int k = 0; k < 3; k++) {
                    term4 += Gamma3[m][i][k] * Gamma3[k][j][m]; 
                }
            }
            R3[i][j] = term1 - term2 + term3 - term4;
        }
    }

}

