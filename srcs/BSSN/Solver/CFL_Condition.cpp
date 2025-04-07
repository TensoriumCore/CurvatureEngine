#include <Geodesics.h>
#include <algorithm>

/**
 * Compute the maximum speed in the grid
 * @return the maximum speed in the grid
 * need to be improved
 */

float Grid::computeMaxSpeed() {
    float maxSpeed = 0.0;
    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int k = 1; k < NZ - 1; k++) {
                Cell2D &cell = globalGrid[i][j][k];
                float betaNorm = std::sqrt(cell.gauge.beta[0]*cell.gauge.beta[0] +
                                            cell.gauge.beta[1]*cell.gauge.beta[1] +
                                            cell.gauge.beta[2]*cell.gauge.beta[2]);
                float localSpeed = std::fabs(cell.gauge.alpha) + betaNorm;
                if (localSpeed > maxSpeed) {
                    maxSpeed = localSpeed;
                }
            }
        }
    }
    return maxSpeed;
}

float Grid::computeCFL_dt(float CFL) {
    float dx_min = std::min({DX, DY, DZ});
    float maxSpeed = computeMaxSpeed();
    
    if (maxSpeed < 1e-10) {
        return 1e-10;
    }
    
    return CFL * dx_min / maxSpeed;
}
