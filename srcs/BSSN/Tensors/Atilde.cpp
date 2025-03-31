#include <Geodesics.h>

void GridTensor::compute_Atilde(Grid &grid_obj, int i, int j, int k) {
    Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
    double Ktrace = 0.0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            Ktrace += cell.geom.gamma_inv[a][b] * cell.curv.K[a][b];
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double Aij = cell.curv.K[a][b] - (1.0/3.0) * Ktrace * cell.geom.gamma[a][b];
            cell.atilde.Atilde[a][b] = std::sqrt(cell.chi) * Aij; 
        }
    }
}

