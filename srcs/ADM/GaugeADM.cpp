#include <Geodesics.h>

void Grid::compute_gauge_derivatives(int i, int j, int k, double &d_alpha_dt, double d_beta_dt[3]) {
    Grid::Cell2D &cell = globalGrid[i][j][k];
    double gammaLocal[3][3], KLocal[3][3];
	GridTensor gridTensor;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaLocal[a][b] = cell.gamma[a][b];
            KLocal[a][b]     = cell.K[a][b];
        }
    }

    double gammaInv[3][3];
    bool ok = invert_3x3(gammaLocal, gammaInv);
    if (!ok) {
        d_alpha_dt = 0.0;
        for (int m = 0; m < 3; m++) d_beta_dt[m] = 0.0;
        return;
	}

    double Ktrace = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * KLocal[a][b];
		}
    }

    double lambda = 1.0 / (1.0 + 2.0 * Ktrace * Ktrace);
    d_alpha_dt = -2.0 * cell.alpha * Ktrace * lambda;

    double eta = 2.0 / (1.0 + std::fabs(Ktrace));
    double d_Gamma_dt[3] = {0.0, 0.0, 0.0}; 

    double tildeGamma[3];
    gridTensor.compute_tildeGamma(i, j, k, tildeGamma); 
    gridTensor.compute_dt_tildeGamma(i, j, k, d_Gamma_dt);

    for (int m = 0; m < 3; m++) {
        d_beta_dt[m] = 3.0 / 4.0 * d_Gamma_dt[m] - eta * cell.beta[m];
    }
}

std::vector<std::vector<std::vector<double>>> hamiltonianGrid;
void Grid::initialize_grid() {
    globalGrid.resize(NX, std::vector<std::vector<Cell2D>>(NY, std::vector<Cell2D>(NZ)));
    hamiltonianGrid.resize(NX, std::vector<std::vector<double>>(NY, std::vector<double>(NZ, 0.0)));
}


void Grid::compute_constraints(Grid &grid_obj, int i, int j, int k, double &hamiltonian, double momentum[3]) {
    double Ricci[3][3];
    compute_ricci_3D(grid_obj, i, j, k, Ricci);
    Cell2D &cell = globalGrid[i][j][k];
    double gammaLocal[3][3], gammaInv[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaLocal[a][b] = cell.gamma[a][b];
        }
    }
    invert_3x3(gammaLocal, gammaInv);
    double R = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            R += gammaInv[a][b] * Ricci[a][b];
        }
    }
    double Ktrace = 0.0;
    double KK = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * cell.K[a][b];
            for (int c = 0; c < 3; c++) {
                KK += cell.K[a][b] * gammaInv[a][c] * cell.K[c][b];
            }
        }
    }
    hamiltonian = R + Ktrace * Ktrace - KK;
    hamiltonianGrid[i][j][k] = hamiltonian;
    double Kup[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double s = 0.0;
            for (int c = 0; c < 3; c++) {
                s += gammaInv[a][c] * cell.K[c][b];
            }
            Kup[a][b] = s;
        }
    }
    for (int i_comp = 0; i_comp < 3; i_comp++) {
        double dCoord = (i_comp == 0 ? DX : (i_comp == 1 ? DY : DZ));
        momentum[i_comp] = 0.0;
        {
            int ip = i, im = i, jp = j, jm = j, kp = k, km = k;
            if (i_comp == 0) { ip = i + 1; im = i - 1; }
            if (i_comp == 1) { jp = j + 1; jm = j - 1; }
            if (i_comp == 2) { kp = k + 1; km = k - 1; }
            if (ip >= 1 && ip < NX && im >= 0 && jm >= 0 && jp < NY && km >= 0 && kp < NZ) {
                double KupPlus[3][3], KupMinus[3][3];
                {
                    Cell2D &cP = globalGrid[ip][jp][kp];
                    double gP[3][3], gInvP[3][3];
                    for (int aa = 0; aa < 3; aa++) {
                        for (int bb = 0; bb < 3; bb++) {
                            gP[aa][bb] = cP.gamma[aa][bb];
                        }
                    }
                    invert_3x3(gP, gInvP);
                    for (int aa = 0; aa < 3; aa++) {
                        for (int bb = 0; bb < 3; bb++) {
                            double tmp = 0.0;
                            for (int cc = 0; cc < 3; cc++) {
                                tmp += gInvP[aa][cc] * cP.K[cc][bb];
                            }
                            KupPlus[aa][bb] = tmp;
                        }
                    }
                }
                {
                    Cell2D &cM = globalGrid[im][jm][km];
                    double gM[3][3], gInvM[3][3];
                    for (int aa = 0; aa < 3; aa++) {
                        for (int bb = 0; bb < 3; bb++) {
                            gM[aa][bb] = cM.gamma[aa][bb];
                        }
                    }
                    invert_3x3(gM, gInvM);
                    for (int aa = 0; aa < 3; aa++) {
                        for (int bb = 0; bb < 3; bb++) {
                            double tmp = 0.0;
                            for (int cc = 0; cc < 3; cc++) {
                                tmp += gInvM[aa][cc] * cM.K[cc][bb];
                            }
                            KupMinus[aa][bb] = tmp;
                        }
                    }
                }
                double partialUp = (KupPlus[i_comp][i_comp] - KupMinus[i_comp][i_comp]) / (2.0 * dCoord);
                momentum[i_comp] += partialUp;
            }
        }
        {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (int j_ = 0; j_ < 3; j_++) {
                for (int m = 0; m < 3; m++) {
                    sum1 += cell.Christoffel[j_][j_][m] * Kup[m][i_comp];
                    sum2 += cell.Christoffel[m][j_][i_comp] * Kup[j_][m];
                }
            }
            momentum[i_comp] += (sum1 - sum2);
        }
        {
            int ip = i, im = i, jp = j, jm = j, kp = k, km = k;
            if (i_comp == 0) { ip = i + 1; im = i - 1; }
            if (i_comp == 1) { jp = j + 1; jm = j - 1; }
            if (i_comp == 2) { kp = k + 1; km = k - 1; }
            if (ip >= 1 && ip < NX && im >= 0 && jm >= 0 && jp < NY && km >= 0 && kp < NZ) {
                double KtraceP = 0.0, KtraceM = 0.0;
                {
                    Cell2D &cP = globalGrid[ip][jp][kp];
                    double gP[3][3], gInvP[3][3];
                    for (int aa = 0; aa < 3; aa++) {
                        for (int bb = 0; bb < 3; bb++) {
                            gP[aa][bb] = cP.gamma[aa][bb];
                        }
                    }
                    invert_3x3(gP, gInvP);
                    for (int aa = 0; aa < 3; aa++) {
                        for (int bb = 0; bb < 3; bb++) {
                            KtraceP += gInvP[aa][bb] * cP.K[aa][bb];
                        }
                    }
                }
                {
                    Cell2D &cM = globalGrid[im][jm][km];
                    double gM[3][3], gInvM[3][3];
                    for (int aa = 0; aa < 3; aa++) {
                        for (int bb = 0; bb < 3; bb++) {
                            gM[aa][bb] = cM.gamma[aa][bb];
                        }
                    }
                    invert_3x3(gM, gInvM);
                    for (int aa = 0; aa < 3; aa++) {
                        for (int bb = 0; bb < 3; bb++) {
                            KtraceM += gInvM[aa][bb] * cM.K[aa][bb];
                        }
                    }
                }
                double dKtrace = (KtraceP - KtraceM) / (2.0 * dCoord);
                momentum[i_comp] -= dKtrace;
            }
        }
    }
}
