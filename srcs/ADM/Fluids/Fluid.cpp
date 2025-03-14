#include <Geodesics.h>

/* here is the implementation of the fluid update functions 
 * in the Grid class
 *
 * This is a simple implematations of a perfect fluid
 * with a constant speed of sound
 * */

void Grid::compute_fluid_derivatives(int i, int j, int k) {
    Cell2D &cell = globalGrid[i][j][k];

    double dpdx = (globalGrid[i+1][j][k].p - globalGrid[i-1][j][k].p) / (2.0 * DX);
    double dpdy = (globalGrid[i][j+1][k].p - globalGrid[i][j-1][k].p) / (2.0 * DY);
    double dpdz = (globalGrid[i][j][k+1].p - globalGrid[i][j][k-1].p) / (2.0 * DZ);

    cell.vx += -dpdx * DT / cell.rho;
    cell.vy += -dpdy * DT / cell.rho;
    cell.vz += -dpdz * DT / cell.rho;

    double drho_dt = -(cell.vx * dpdx + cell.vy * dpdy + cell.vz * dpdz);
    cell.rho += drho_dt * DT;
}


void Grid::update_fluid_velocity(int i, int j, int k, double dt) {
    Cell2D &cell = globalGrid[i][j][k];

    if (cell.rho < 1e-10) return;

    double pressure_gradient_x = (globalGrid[i+1][j][k].p - globalGrid[i-1][j][k].p) / (2.0 * DX);
    double pressure_gradient_y = (globalGrid[i][j+1][k].p - globalGrid[i][j-1][k].p) / (2.0 * DY);
    double pressure_gradient_z = (globalGrid[i][j][k+1].p - globalGrid[i][j][k-1].p) / (2.0 * DZ);

    double christoffel_x = 0.0, christoffel_y = 0.0, christoffel_z = 0.0;
    for (int b = 0; b < 3; b++) {
        for (int c = 0; c < 3; c++) {
            christoffel_x += cell.Christoffel[0][b][c] * (b == 0 ? cell.vx : (b == 1 ? cell.vy : cell.vz)) *
                                                         (c == 0 ? cell.vx : (c == 1 ? cell.vy : cell.vz));
            christoffel_y += cell.Christoffel[1][b][c] * (b == 0 ? cell.vx : (b == 1 ? cell.vy : cell.vz)) *
                                                         (c == 0 ? cell.vx : (c == 1 ? cell.vy : cell.vz));
            christoffel_z += cell.Christoffel[2][b][c] * (b == 0 ? cell.vx : (b == 1 ? cell.vy : cell.vz)) *
                                                         (c == 0 ? cell.vx : (c == 1 ? cell.vy : cell.vz));
        }
    }

    cell.vx += -dt * (pressure_gradient_x / cell.rho + christoffel_x);
    cell.vy += -dt * (pressure_gradient_y / cell.rho + christoffel_y);
    cell.vz += -dt * (pressure_gradient_z / cell.rho + christoffel_z);
	/* printf("vx = %f, vy = %f, vz = %f\n", cell.vx, cell.vy, cell.vz); */
}
