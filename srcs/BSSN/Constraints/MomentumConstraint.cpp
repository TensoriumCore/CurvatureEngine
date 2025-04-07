#include <Geodesics.h>

std::vector<std::vector<std::vector<float>>> hamiltonianGrid;
void Grid::initialize_grid() {
    globalGrid.resize(NX, std::vector<std::vector<Cell2D>>(NY, std::vector<Cell2D>(NZ)));
    hamiltonianGrid.resize(NX, std::vector<std::vector<float>>(NY, std::vector<float>(NZ, 0.0)));
}


void Grid::compute_constraints(Grid &grid_obj, int i, int j, int k, float &hamiltonian, float momentum[3]) {
    Cell2D &cell = globalGrid[i][j][k];
    float R = 0.0;
	GridTensor grid_tensor_obj;
	R = compute_ricci_scalar(grid_obj, i, j, k);
    float Ktrace = 0.0;
    float KK = 0.0;
#pragma omp simd collapse(4)
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					KK += cell.curv.K[i][j] * 
						cell.geom.gamma_inv[i][k] * 
						cell.geom.gamma_inv[j][l] * 
						cell.curv.K[k][l];
				}
			}
		}
	}

	Log log_obj;
    hamiltonian = R + Ktrace * Ktrace - KK;
    hamiltonianGrid[i][j][k] = hamiltonian;
	cell.matter.hamiltonian = hamiltonian;
	for(int i_comp=0; i_comp<3; i_comp++){
		cell.matter.momentum[i_comp] = grid_tensor_obj.compute_momentum_i(grid_obj, i, j, k, i_comp);
	}
}
