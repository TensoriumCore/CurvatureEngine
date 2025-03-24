#include <Geodesics.h>

/* void Grid::update_energy_momentum_tensor(int i, int j, int k) { */
/*     Cell2D &cell = globalGrid[i][j][k]; */
/*  */
/*     double v2 = cell.vx * cell.vx + cell.vy * cell.vy + cell.vz * cell.vz; */
/*     double gamma_lorentz = 1.0 / sqrt(1.0 - v2);  */
/*  */
/*     double rho_h = cell.rho + cell.p;  */
/*  */
/*     for (int mu = 0; mu < 4; mu++) { */
/*         for (int nu = 0; nu < 4; nu++) { */
/*             cell.T[mu][nu] = 0.0; */
/*         } */
/*     } */
/*  */
/*     cell.T[0][0] = rho_h * gamma_lorentz * gamma_lorentz; */
/*     cell.T[0][1] = rho_h * gamma_lorentz * gamma_lorentz * cell.vx; */
/*     cell.T[0][2] = rho_h * gamma_lorentz * gamma_lorentz * cell.vy; */
/*     cell.T[0][3] = rho_h * gamma_lorentz * gamma_lorentz * cell.vz; */
/*  */
/*     cell.T[1][0] = cell.T[0][1]; */
/*     cell.T[2][0] = cell.T[0][2]; */
/*     cell.T[3][0] = cell.T[0][3]; */
/*  */
/*     cell.T[1][1] = rho_h * gamma_lorentz * gamma_lorentz * cell.vx * cell.vx + cell.p; */
/*     cell.T[1][2] = rho_h * gamma_lorentz * gamma_lorentz * cell.vx * cell.vy; */
/*     cell.T[1][3] = rho_h * gamma_lorentz * gamma_lorentz * cell.vx * cell.vz; */
/*  */
/*     cell.T[2][1] = cell.T[1][2]; */
/*     cell.T[2][2] = rho_h * gamma_lorentz * gamma_lorentz * cell.vy * cell.vy + cell.p; */
/*     cell.T[2][3] = rho_h * gamma_lorentz * gamma_lorentz * cell.vy * cell.vz; */
/*  */
/*     cell.T[3][1] = cell.T[1][3]; */
/*     cell.T[3][2] = cell.T[2][3]; */
/*     cell.T[3][3] = rho_h * gamma_lorentz * gamma_lorentz * cell.vz * cell.vz + cell.p; */
/* } */
/*  */
/* void Grid::compute_energy_momentum_evolution(int i, int j, int k, double dt) { */
/*     Cell2D &cell = globalGrid[i][j][k]; */
/*  */
/*     double dTdx[4][4], dTdy[4][4], dTdz[4][4]; */
/*     for (int mu = 0; mu < 4; mu++) { */
/*         for (int nu = 0; nu < 4; nu++) { */
/*             dTdx[mu][nu] = (globalGrid[i+1][j][k].T[mu][nu] - globalGrid[i-1][j][k].T[mu][nu]) / (2.0 * DX); */
/*             dTdy[mu][nu] = (globalGrid[i][j+1][k].T[mu][nu] - globalGrid[i][j-1][k].T[mu][nu]) / (2.0 * DY); */
/*             dTdz[mu][nu] = (globalGrid[i][j][k+1].T[mu][nu] - globalGrid[i][j][k-1].T[mu][nu]) / (2.0 * DZ); */
/*         } */
/*     } */
/*  */
/*     double ChristoffelTerm[4][4] = {0.0}; */
/*     for (int nu = 0; nu < 4; nu++) { */
/*         for (int mu = 0; mu < 4; mu++) { */
/*             for (int alpha = 0; alpha < 4; alpha++) { */
/*                 for (int beta = 0; beta < 4; beta++) { */
/*                     ChristoffelTerm[mu][nu] += cell.Christoffel[nu][alpha][beta] * cell.T[alpha][beta]; */
/*                 } */
/*             } */
/*         } */
/*     } */
/*  */
/*     double dtT[4][4]; */
/*     for (int mu = 0; mu < 4; mu++) { */
/*         for (int nu = 0; nu < 4; nu++) { */
/*             double divergence = dTdx[mu][nu] + dTdy[mu][nu] + dTdz[mu][nu]; */
/*             dtT[mu][nu] = -divergence - ChristoffelTerm[mu][nu]; */
/*         } */
/*     } */
/*  */
/*     for (int mu = 0; mu < 4; mu++) { */
/*         for (int nu = 0; nu < 4; nu++) { */
/*             cell.T[mu][nu] += dt * dtT[mu][nu]; */
/*         } */
/*     } */
/* } */
/*  */
