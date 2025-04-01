/* #include <Geodesics.h> */
/*  */
/* double Grid::compute_ADM_mass() { */
/*     double mass = 0.0; */
/*  */
/*     int i_faces[] = {1, NX - 2}; */
/*     int j_faces[] = {1, NY - 2}; */
/*     int k_faces[] = {1, NZ - 2}; */
/*  */
/*     double inv_16pi = 1.0 / (16.0 * M_PI); */
/*  */
/* 	printf("gamma_xx at (NX-2, NY/2, NZ/2): %e\n", globalGrid[NX-2][NY/2][NZ/2].geom.gamma[0][0]); */
/*     for (int f = 0; f < 2; f++) { */
/*         int i = i_faces[f]; */
/*         int sign = (f == 0) ? -1 : 1; */
/*         #pragma omp parallel for reduction(+:mass) collapse(2) */
/*         for (int j = 1; j < NY - 1; j++) { */
/*             for (int k = 1; k < NZ - 1; k++) { */
/*                 Cell2D &cell     = globalGrid[i][j][k]; */
/*                 Cell2D &cell_adj = globalGrid[i + sign][j][k]; */
/*  */
/*                 double dgamma_xx_dx = (cell.geom.gamma[0][0] - cell_adj.geom.gamma[0][0]) / DX; */
/*                 double dgamma_yy_dx = (cell.geom.gamma[1][1] - cell_adj.geom.gamma[1][1]) / DX; */
/*                 double dgamma_zz_dx = (cell.geom.gamma[2][2] - cell_adj.geom.gamma[2][2]) / DX; */
/*  */
/*                 double flux = dgamma_xx_dx - (dgamma_yy_dx + dgamma_zz_dx); */
/*                 mass += flux * DY * DZ; */
/*             } */
/*         } */
/*     } */
/*  */
/*     for (int f = 0; f < 2; f++) { */
/*         int j = j_faces[f]; */
/*         int sign = (f == 0) ? -1 : 1; */
/*         #pragma omp parallel for reduction(+:mass) collapse(2) */
/*         for (int i = 1; i < NX - 1; i++) { */
/*             for (int k = 1; k < NZ - 1; k++) { */
/*                 Cell2D &cell     = globalGrid[i][j][k]; */
/*                 Cell2D &cell_adj = globalGrid[i][j + sign][k]; */
/*  */
/*                 double dgamma_yy_dy = (cell.geom.gamma[1][1] - cell_adj.geom.gamma[1][1]) / DY; */
/*                 double dgamma_xx_dy = (cell.geom.gamma[0][0] - cell_adj.geom.gamma[0][0]) / DY; */
/*                 double dgamma_zz_dy = (cell.geom.gamma[2][2] - cell_adj.geom.gamma[2][2]) / DY; */
/*  */
/*                 double flux = dgamma_yy_dy - (dgamma_xx_dy + dgamma_zz_dy); */
/*                 mass += flux * DX * DZ; */
/*             } */
/*         } */
/*     } */
/*  */
/*     for (int f = 0; f < 2; f++) { */
/*         int k = k_faces[f]; */
/*         int sign = (f == 0) ? -1 : 1; */
/*         #pragma omp parallel for reduction(+:mass) collapse(2) */
/*         for (int i = 1; i < NX - 1; i++) { */
/*             for (int j = 1; j < NY - 1; j++) { */
/*                 Cell2D &cell     = globalGrid[i][j][k]; */
/*                 Cell2D &cell_adj = globalGrid[i][j][k + sign]; */
/*  */
/*                 double dgamma_zz_dz = (cell.geom.gamma[2][2] - cell_adj.geom.gamma[2][2]) / DZ; */
/*                 double dgamma_xx_dz = (cell.geom.gamma[0][0] - cell_adj.geom.gamma[0][0]) / DZ; */
/*                 double dgamma_yy_dz = (cell.geom.gamma[1][1] - cell_adj.geom.gamma[1][1]) / DZ; */
/*  */
/*                 double flux = dgamma_zz_dz - (dgamma_xx_dz + dgamma_yy_dz); */
/*                 mass += flux * DX * DY; */
/*             } */
/*         } */
/*     } */
/*  */
/*     mass *= inv_16pi; */
/*  */
/*     printf("Approx. ADM mass (6 faces) = %.6e\n", mass); */
/*     return mass; */
/* } */
