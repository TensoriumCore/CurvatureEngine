#include <Geodesics.h>
#include <cassert>


void Grid::compute_time_derivatives(Grid &grid_obj, int i, int j, int k)
{
    Cell2D &cell = globalGrid[i][j][k];
    double alpha = cell.alpha;
    double beta[3] = { cell.beta[0], cell.beta[1], cell.beta[2] };
    double gammaLocal[3][3];
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            gammaLocal[a][b] = cell.gamma[a][b];
            cell.K[a][b]     = cell.K[a][b];
        }
    }

    /* double cell.gamma_inv[3][3]; */
    /* invert_3x3(gammaLocal, cell.gamma_inv); */
    /*  */
    double Ktrace = 0.0;
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            Ktrace += cell.gamma_inv[a][b] * cell.K[a][b];
        }
    }

    GridTensor gridTensor;
    double Gamma[3][3][3];
    gridTensor.compute_christoffel_3D(grid_obj, i, j, k, Gamma);

    double Ricci[3][3];
    compute_ricci_3D(grid_obj, i, j, k, Ricci);

    
	auto partial_m_alpha = [&](int ii, int jj, int kk, int m) {
		auto alpha_accessor = [](const Grid::Cell2D &cell) {
			return cell.alpha;
		};

		if (m == 0) return partialX(grid_obj, ii, jj, kk, alpha_accessor);
		if (m == 1) return partialY(grid_obj, ii, jj, kk, alpha_accessor);
		if (m == 2) return partialZ(grid_obj, ii, jj, kk, alpha_accessor);
		return 0.0;
	};

    auto D_iD_j_alpha = [&](int a, int b) {
        double secondPart = second_partial_alpha(grid_obj, i, j, k, a, b);
        double sumG = 0.0;
        for(int m = 0; m < 3; m++){
            sumG += Gamma[m][a][b] * partial_m_alpha(i, j, k, m);
        }
        return (secondPart - sumG);
    };

	auto partial_m_K = [&](int a, int b, int dim) {
		return partial_m(grid_obj, i, j, k, dim, [&](const Grid::Cell2D &cell) {
				return cell.K[a][b];
				});
	};

	auto partial_m_gamma = [&](int a, int b, int dim) {
		return partial_m(grid_obj, i, j, k, dim, [&](const Grid::Cell2D &cell) {
				return cell.gamma[a][b];
				});
	};
    
	auto partial_j_beta_comp = [&](int dimSpace, int compBeta) {
		return partial_m(grid_obj, i, j, k, dimSpace, [&](const Grid::Cell2D &cell) {
				return cell.beta[compBeta];
				});
	};

    double partialBeta[3][3];
    for(int jDim = 0; jDim < 3; jDim++){
        for(int iComp = 0; iComp < 3; iComp++){
            partialBeta[jDim][iComp] = partial_j_beta_comp(jDim, iComp);

        }
    }
    double dtGamma[3][3];
    for(int a = 0; a < 3; a++){
        for(int b = 0; b < 3; b++){
            double rhs = -2.0 * alpha * cell.K[a][b];

            double advection = 0.0;
            for(int m = 0; m < 3; m++){
                advection += beta[m] * partial_m_gamma(a, b, m);
            }

            double shiftTerm = 0.0;
            for(int m = 0; m < 3; m++){
                shiftTerm += gammaLocal[a][m] * partialBeta[b][m];
                shiftTerm += gammaLocal[b][m] * partialBeta[a][m];
            }

            dtGamma[a][b] = rhs + advection + shiftTerm;
        }
    }

if (i == NX/2 && j == NY/2 && k == NZ/2) {
    printf("K_00 = %e, K_11 = %e, K_22 = %e\n", cell.K[0][0], cell.K[1][1], cell.K[2][2]);
}

    double dtK[3][3];
    for(int a = 0; a < 3; a++){
        for(int b = 0; b < 3; b++){
            double dDaDb_alpha = - D_iD_j_alpha(a, b);

            double KContra[3][3];
            for(int m = 0; m < 3; m++){
                for(int n = 0; n < 3; n++){
                    double sum_ = 0.0;
                    for(int r = 0; r < 3; r++){
                        for(int s = 0; s < 3; s++){
                            sum_ += cell.gamma_inv[m][r] * cell.gamma_inv[n][s] * cell.K[r][s];
                        }
                    }
                    KContra[m][n] = sum_;
                }
            }

            double KtraceLocal = 0.0;
            for(int x = 0; x < 3; x++){
                for(int y = 0; y < 3; y++){
                    KtraceLocal += cell.gamma_inv[x][y] * cell.K[x][y];
                }
            }

            double RicciTerm  = Ricci[a][b];
            double KKab       = KtraceLocal * cell.K[a][b];

            double KaM_KM_b = 0.0;
            for(int m = 0; m < 3; m++){
                KaM_KM_b += cell.K[a][m] * KContra[m][b];
            }

            double termB = alpha * (RicciTerm + KKab - 2.0 * KaM_KM_b);

            double advK = 0.0;
            for(int m = 0; m < 3; m++){
                advK += beta[m] * partial_m_K(a, b, m);
            }

            double shiftK = 0.0;
            for(int m = 0; m < 3; m++){
                shiftK += cell.K[m][b] * partialBeta[a][m];
                shiftK += cell.K[a][m] * partialBeta[b][m];
            }

            double lieBetaK = advK + shiftK;

            dtK[a][b] = dDaDb_alpha + termB + lieBetaK;
        }
    }

    for(int a = 0; a < 3; a++){
        for(int b = 0; b < 3; b++){
            cell.dgt[a][b] = dtGamma[a][b]; 
            cell.dKt[a][b] = dtK[a][b];  
        }
    }

    double norm_dtK = 0.0;
    for(int a = 0; a < 3; a++){
        for(int b = 0; b < 3; b++){
            norm_dtK += dtK[a][b] * dtK[a][b];
        }
    }
    norm_dtK = sqrt(norm_dtK / 9.0);
    /* printf("norm(dtK) = %f\n", norm_dtK); */
}
