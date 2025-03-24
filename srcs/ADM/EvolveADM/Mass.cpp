#include <Geodesics.h>

double Grid::compute_ADM_mass() {
    double mass = 0.0;

    int iSurf = NX - 2;

    #pragma omp parallel for reduction(+:mass)
    for (int j = 1; j < NY - 1; j++) {
        for (int k = 1; k < NZ - 1; k++) {

            Cell2D &cell_i  = globalGrid[iSurf][j][k];
            Cell2D &cell_im = globalGrid[iSurf - 1][j][k];

            double dgamma_xx_dx = (cell_i.geom.gamma[0][0] - cell_im.geom.gamma[0][0]) / DX;
            double dgamma_yy_dx = (cell_i.geom.gamma[1][1] - cell_im.geom.gamma[1][1]) / DX;
            double dgamma_zz_dx = (cell_i.geom.gamma[2][2] - cell_im.geom.gamma[2][2]) / DX;

            double flux = dgamma_xx_dx - (dgamma_yy_dx + dgamma_zz_dx);

            double dS = DY * DZ;

            mass += flux * dS;
        }
    }

    mass *= 1.0 / (16.0 * M_PI);

    printf("Approx. ADM mass (face x=NX-2) = %.6e\n", mass);
    return mass;
}
