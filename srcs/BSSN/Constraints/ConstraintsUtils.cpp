#include <Geodesics.h>

float Grid::KUpAt(Grid &grid, int ip, int jp, int kp, int j_up, int i_low)
{
    if(ip<0 || ip>=NX || jp<0 || jp>=NY || kp<0 || kp>=NZ) {
        return 0.0;
    }

    Cell2D &cell = grid.globalGrid[ip][jp][kp];
    float gtmp[3][3];
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            gtmp[a][b] = cell.geom.gamma[a][b];
        }
    }

    float val = 0.0;
    for(int c=0; c<3; c++){
        val += cell.geom.gamma_inv[j_up][c] * cell.curv.K[c][i_low];
    }
    return val;
}

float Grid::partialX_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low)
{
    float valP = KUpAt(grid, i+1, j, k, j_up, i_low);
    float valM = KUpAt(grid, i-1, j, k, j_up, i_low);
    return (valP - valM) / (2.0 * DX);
}

float Grid::partialY_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low)
{
    float valP = KUpAt(grid, i, j+1, k, j_up, i_low);
    float valM = KUpAt(grid, i, j-1, k, j_up, i_low);
    return (valP - valM) / (2.0 * DY);
}

float Grid::partialZ_KUp(Grid &grid, int i, int j, int k, int j_up, int i_low)
{
    float valP = KUpAt(grid, i, j, k+1, j_up, i_low);
    float valM = KUpAt(grid, i, j, k-1, j_up, i_low);
    return (valP - valM) / (2.0 * DZ);
}

float divKUp_i(Grid &grid, int i, int j, int k, int i_comp)
{
    float px = grid.partialX_KUp(grid, i, j, k, 0, i_comp);
    float py = grid.partialY_KUp(grid, i, j, k, 1, i_comp);
    float pz = grid.partialZ_KUp(grid, i, j, k, 2, i_comp);
    return px + py + pz;
}

float GridTensor::christoffelTerm(Grid &grid, int i, int j, int k, int i_comp)
{
    float KupLocal[3][3];
    {
        float tmpG[3][3];
#pragma omp simd
        for(int a=0; a<3; a++){
            for(int b=0; b<3; b++){
                tmpG[a][b] = grid.getCell(i,j,k).geom.gamma[a][b];
            }
        }
        for(int aa=0; aa<3; aa++){
            for(int bb=0; bb<3; bb++){
                float sum=0.0;
                for(int cc=0; cc<3; cc++){
					sum += grid.getCell(i,j,k).geom.gamma_inv[aa][cc] * grid.getCell(i,j,k).curv.K[cc][bb];
                }
                KupLocal[aa][bb] = sum;
            }
        }
    }
    float sum1=0.0, sum2=0.0;
    for(int j_=0; j_<3; j_++){
        for(int m=0; m<3; m++){
			sum1 += grid.getCell(i,j,k).conn.Christoffel[m][j_][j_] * KupLocal[m][i_comp];
			sum2 += grid.getCell(i,j,k).conn.Christoffel[m][j_][i_comp] * KupLocal[j_][m];
        }
    }
    return sum1 - sum2;
}


float Grid::computeTraceK(Grid &grid, int i, int j, int k)
{
    Cell2D &cell = grid.globalGrid[i][j][k];
    float gTmp[3][3];
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            gTmp[a][b] = cell.geom.gamma[a][b];
        }
    }
    float trace = 0.0;
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            trace += cell.geom.gamma_inv[a][b]*cell.curv.K[a][b];
        }
    }
    return trace;
}

float partialX_Ktrace(Grid &grid, int i, int j, int k)
{
    float Kp = grid.computeTraceK(grid, i+1,j,k);
    float Km = grid.computeTraceK(grid, i-1,j,k);
    return (Kp - Km)/(2.0*DX);
}

float partialY_Ktrace(Grid &grid, int i, int j, int k)
{
	float Kp = grid.computeTraceK(grid, i,j+1,k);
	float Km = grid.computeTraceK(grid, i,j-1,k);
	return (Kp - Km)/(2.0*DY);
}

float partialZ_Ktrace(Grid &grid, int i, int j, int k)
{
	float Kp = grid.computeTraceK(grid, i,j,k+1);
	float Km = grid.computeTraceK(grid, i,j,k-1);
	return (Kp - Km)/(2.0*DZ);
}

float partial_i_Ktrace(Grid &grid, int i, int j, int k, int i_comp)
{
    if(i_comp==0) return partialX_Ktrace(grid, i,j,k);
    if(i_comp==1) return partialY_Ktrace(grid, i,j,k);
    if(i_comp==2) return partialZ_Ktrace(grid, i,j,k);
    return 0.0;
}

float GridTensor::compute_momentum_i(Grid &grid, int i, int j, int k, int i_comp)
{
    float divKUpVal = divKUp_i(grid, i, j, k, i_comp);
    float christofVal = christoffelTerm(grid, i, j, k, i_comp);
    float dK = partial_i_Ktrace(grid, i, j, k, i_comp);
    return divKUpVal + christofVal - dK;
}

float Grid::compute_ricci_scalar(Grid &grid, int i, int j, int k)
{
	Cell2D &cell = grid.globalGrid[i][j][k];
	Log log_obj;
	float gTmp[3][3];
	for(int a=0; a<3; a++){
		for(int b=0; b<3; b++){
			gTmp[a][b] = cell.geom.gamma[a][b];
		}
	}
	float R = 0.0;
	for(int a=0; a<3; a++){
		for(int b=0; b<3; b++){
			R += cell.geom.gamma_inv[a][b]*cell.geom.Ricci[a][b];
		}
	}
	return R;
} 


