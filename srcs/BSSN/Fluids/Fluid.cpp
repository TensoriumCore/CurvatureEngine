#include <Geodesics.h>

/* here is the implementation of the fluid update functions 
 * in the Grid class
 *
 * This is a simple implematations of a perfect fluid
 * with a constant speed of sound
 * */

void Grid::compute_fluid_derivatives(int i, int j, int k) {
    Cell2D &cell = globalGrid[i][j][k];

    double dpdx = (globalGrid[i+1][j][k].matter.p - globalGrid[i-1][j][k].matter.p) / (2.0 * DX);
    double dpdy = (globalGrid[i][j+1][k].matter.p - globalGrid[i][j-1][k].matter.p) / (2.0 * DY);
    double dpdz = (globalGrid[i][j][k+1].matter.p - globalGrid[i][j][k-1].matter.p) / (2.0 * DZ);

    cell.matter.vx += -dpdx * DT / cell.matter.rho;
    cell.matter.vy += -dpdy * DT / cell.matter.rho;
    cell.matter.vz += -dpdz * DT / cell.matter.rho;

    double drho_dt = -(cell.matter.vx * dpdx + cell.matter.vy * dpdy + cell.matter.vz * dpdz);
    cell.matter.rho += drho_dt * DT;
}


void Grid::update_fluid_velocity(int i, int j, int k, double dt) {
    Cell2D &cell = globalGrid[i][j][k];

    if (cell.matter.rho < 1e-10) return;

    double pressure_gradient_x = (globalGrid[i+1][j][k].matter.p - globalGrid[i-1][j][k].matter.p) / (2.0 * DX);
    double pressure_gradient_y = (globalGrid[i][j+1][k].matter.p - globalGrid[i][j-1][k].matter.p) / (2.0 * DY);
    double pressure_gradient_z = (globalGrid[i][j][k+1].matter.p - globalGrid[i][j][k-1].matter.p) / (2.0 * DZ);

    double christoffel_x = 0.0, christoffel_y = 0.0, christoffel_z = 0.0;
    for (int b = 0; b < 3; b++) {
        for (int c = 0; c < 3; c++) {
            christoffel_x += cell.conn.Christoffel[0][b][c] * (b == 0 ? cell.matter.vx : (b == 1 ? cell.matter.vy : cell.matter.vz)) *
                                                         (c == 0 ? cell.matter.vx : (c == 1 ? cell.matter.vy : cell.matter.vz));
            christoffel_y += cell.conn.Christoffel[1][b][c] * (b == 0 ? cell.matter.vx : (b == 1 ? cell.matter.vy : cell.matter.vz)) *
                                                         (c == 0 ? cell.matter.vx : (c == 1 ? cell.matter.vy : cell.matter.vz));
            christoffel_z += cell.conn.Christoffel[2][b][c] * (b == 0 ? cell.matter.vx : (b == 1 ? cell.matter.vy : cell.matter.vz)) *
                                                         (c == 0 ? cell.matter.vx : (c == 1 ? cell.matter.vy : cell.matter.vz));
        }
    }

    cell.matter.vx += -dt * (pressure_gradient_x / cell.matter.rho + christoffel_x);
    cell.matter.vy += -dt * (pressure_gradient_y / cell.matter.rho + christoffel_y);
    cell.matter.vz += -dt * (pressure_gradient_z / cell.matter.rho + christoffel_z);
	/* printf("matter.vx = %f, vy = %f, vz = %f\n", cell.matter.vx, cell.matter.vy, cell.matter.vz); */
}
