#include <Geodesics.h>
#include <algorithm>

/**
 * Compute the maximum speed in the grid
 * @return the maximum speed in the grid
 * need to be improved
 */

double Grid::computeMaxSpeed() {
    double maxSpeed = 0.0;
    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int k = 1; k < NZ - 1; k++) {
                Cell2D &cell = globalGrid[i][j][k];
                double betaNorm = std::sqrt(cell.beta[0]*cell.beta[0] +
                                            cell.beta[1]*cell.beta[1] +
                                            cell.beta[2]*cell.beta[2]);
                double localSpeed = std::fabs(cell.alpha) + betaNorm;
                if (localSpeed > maxSpeed) {
                    maxSpeed = localSpeed;
                }
            }
        }
    }
    return maxSpeed;
}

double Grid::computeCFL_dt(double CFL) {
    double dx_min = std::min({DX, DY, DZ});
    double maxSpeed = computeMaxSpeed();
    
    if (maxSpeed < 1e-10) {
        return 1e-10;
    }
    
    return CFL * dx_min / maxSpeed;
}
