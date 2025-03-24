#include <Geodesics.h>


std::vector<std::vector<std::vector<double>>> hamiltonianGrid;
void Grid::initialize_grid() {
    globalGrid.resize(NX, std::vector<std::vector<Cell2D>>(NY, std::vector<Cell2D>(NZ)));
    hamiltonianGrid.resize(NX, std::vector<std::vector<double>>(NY, std::vector<double>(NZ, 0.0)));
}

double Grid::KUpAt(Grid &grid, int ip, int jp, int kp, int j_up, int i_low)
{
    if(ip<0 || ip>=NX || jp<0 || jp>=NY || kp<0 || kp>=NZ) {
        return 0.0;
    }

    Cell2D &cell = grid.globalGrid[ip][jp][kp];
    double gtmp[3][3];
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            gtmp[a][b] = cell.geom.gamma[a][b];
        }
    }

    double val = 0.0;
    for(int c=0; c<3; c++){
        val += cell.geom.gamma_inv[j_up][c] * cell.curv.K[c][i_low];
    }
    return val;
}

double Grid::partialX_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low)
{
    double valP = KUpAt(grid, i+1, j, k, j_up, i_low);
    double valM = KUpAt(grid, i-1, j, k, j_up, i_low);
    return (valP - valM) / (2.0 * DX);
}

double Grid::partialY_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low)
{
    double valP = KUpAt(grid, i, j+1, k, j_up, i_low);
    double valM = KUpAt(grid, i, j-1, k, j_up, i_low);
    return (valP - valM) / (2.0 * DY);
}

double Grid::partialZ_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low)
{
    double valP = KUpAt(grid, i, j, k+1, j_up, i_low);
    double valM = KUpAt(grid, i, j, k-1, j_up, i_low);
    return (valP - valM) / (2.0 * DZ);
}

double divKUp_i(Grid &grid, int i, int j, int k, int i_comp)
{
    double px = grid.partialX_KUp(grid, i, j, k, 0, i_comp);
    double py = grid.partialY_KUp(grid, i, j, k, 1, i_comp);
    double pz = grid.partialZ_KUp(grid, i, j, k, 2, i_comp);
    return px + py + pz;
}

double Grid::christoffelTerm(Grid &grid, int i, int j, int k, int i_comp)
{
    Cell2D &cell = grid.globalGrid[i][j][k];
    double KupLocal[3][3];
    {
        double tmpG[3][3];
#pragma omp simd
        for(int a=0; a<3; a++){
            for(int b=0; b<3; b++){
                tmpG[a][b] = cell.geom.gamma[a][b];
            }
        }
        for(int aa=0; aa<3; aa++){
            for(int bb=0; bb<3; bb++){
                double sum=0.0;
                for(int cc=0; cc<3; cc++){
                    sum += cell.geom.gamma_inv[aa][cc]*cell.curv.K[cc][bb];
                }
                KupLocal[aa][bb] = sum;
            }
        }
    }
    double sum1=0.0, sum2=0.0;
    for(int j_=0; j_<3; j_++){
        for(int m=0; m<3; m++){
            sum1 += cell.conn.Christoffel[m][j_][j_] * KupLocal[m][i_comp];
            sum2 += cell.conn.Christoffel[m][j_][i_comp] * KupLocal[j_][m];
        }
    }
    return sum1 - sum2;
}


double Grid::computeTraceK(Grid &grid, int i, int j, int k)
{
    Cell2D &cell = grid.globalGrid[i][j][k];
    double gTmp[3][3];
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            gTmp[a][b] = cell.geom.gamma[a][b];
        }
    }
    double trace = 0.0;
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            trace += cell.geom.gamma_inv[a][b]*cell.curv.K[a][b];
        }
    }
    return trace;
}

double partialX_Ktrace(Grid &grid, int i, int j, int k)
{
    double Kp = grid.computeTraceK(grid, i+1,j,k);
    double Km = grid.computeTraceK(grid, i-1,j,k);
    return (Kp - Km)/(2.0*DX);
}

double partialY_Ktrace(Grid &grid, int i, int j, int k)
{
	double Kp = grid.computeTraceK(grid, i,j+1,k);
	double Km = grid.computeTraceK(grid, i,j-1,k);
	return (Kp - Km)/(2.0*DY);
}

double partialZ_Ktrace(Grid &grid, int i, int j, int k)
{
	double Kp = grid.computeTraceK(grid, i,j,k+1);
	double Km = grid.computeTraceK(grid, i,j,k-1);
	return (Kp - Km)/(2.0*DZ);
}

double partial_i_Ktrace(Grid &grid, int i, int j, int k, int i_comp)
{
    if(i_comp==0) return partialX_Ktrace(grid, i,j,k);
    if(i_comp==1) return partialY_Ktrace(grid, i,j,k);
    if(i_comp==2) return partialZ_Ktrace(grid, i,j,k);
    return 0.0;
}

double compute_momentum_i(Grid &grid, int i, int j, int k, int i_comp)
{
    double divKUpVal = divKUp_i(grid, i, j, k, i_comp);
    double christofVal = grid.christoffelTerm(grid, i, j, k, i_comp);
    double dK = partial_i_Ktrace(grid, i, j, k, i_comp);
    return divKUpVal + christofVal - dK;
}

double Grid::compute_ricci_scalar(Grid &grid, int i, int j, int k)
{
	Cell2D &cell = grid.globalGrid[i][j][k];
	Log log_obj;
	double gTmp[3][3];
	for(int a=0; a<3; a++){
		for(int b=0; b<3; b++){
			gTmp[a][b] = cell.geom.gamma[a][b];
		}
	}
	double R = 0.0;
	for(int a=0; a<3; a++){
		for(int b=0; b<3; b++){
			R += cell.geom.gamma_inv[a][b]*cell.geom.Ricci[a][b];
		}
	}
	return R;
} 


void Grid::compute_constraints(Grid &grid_obj, int i, int j, int k, double &hamiltonian, double momentum[3]) {
    Cell2D &cell = globalGrid[i][j][k];
    double gammaLocal[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaLocal[a][b] = cell.geom.gamma[a][b];
        }
    }
    double R = 0.0;
    /* for (int a = 0; a < 3; a++) { */
    /*     for (int b = 0; b < 3; b++) { */
    /*         R += cell.geom.gamma_inv[a][b] * cell.geom.Ricci[a][b]; */
    /*     } */
    /* } */
	R = compute_ricci_scalar(grid_obj, i, j, k);
    double Ktrace = 0.0;
    double KK = 0.0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					KK += cell.curv.K[i][j] * cell.geom.gamma_inv[i][k] * cell.geom.gamma_inv[j][l] * cell.curv.K[k][l];
				}
			}
		}
	}

	Log log_obj;
    hamiltonian = R + Ktrace * Ktrace - KK;
    hamiltonianGrid[i][j][k] = hamiltonian;
	cell.matter.hamiltonian = hamiltonian;
	for(int i_comp=0; i_comp<3; i_comp++){
		momentum[i_comp] = compute_momentum_i(grid_obj, i, j, k, i_comp);
		cell.matter.momentum[i_comp] = momentum[i_comp];
	}
}
