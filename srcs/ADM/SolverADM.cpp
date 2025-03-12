#include <Geodesics.h>
#include <cassert>



void Grid::compute_time_derivatives(Grid &grid_obj, int i, int j, int k)
{
    Cell2D &cell = globalGrid[i][j][k];
    double alpha = cell.alpha;
    double beta[3] = { cell.beta[0], cell.beta[1], cell.beta[2] };
	GridTensor gridTensor;
    double gammaLocal[3][3], KLocal[3][3];
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            gammaLocal[a][b] = cell.gamma[a][b];
            KLocal[a][b]     = cell.K[a][b];
        }
    }

    double gammaInv[3][3];
    invert_3x3(gammaLocal, gammaInv);

    double Ktrace = 0.0;
    for(int a=0; a<3; a++){
		for(int b=0; b<3; b++){
			Ktrace += gammaInv[a][b]*KLocal[a][b];
		}
	}

    double Gamma[3][3][3];
    gridTensor.compute_christoffel_3D(grid_obj, i,j,k,Gamma);
    
    double Ricci[3][3];
    compute_ricci_3D(grid_obj, i,j,k,Ricci);


	auto partial_m_alpha = [&](int i, int j, int k, int m){
		if(m==0) return partialX_alpha(grid_obj, i, j, k);
		if(m==1) return partialY_alpha(grid_obj, i, j, k);
		if(m==2) return partialZ_alpha(grid_obj, i, j, k);
		return 0.0;
	};

	auto D_iD_j_alpha = [&](int a, int b){
		double secondPart = second_partial_alpha(grid_obj, i,j,k,a,b);
		double sumG = 0.0;
		for(int m=0; m<3; m++){
			sumG += Gamma[m][a][b] * partial_m_alpha(i,j,k,m);
		}
		return secondPart - sumG;
	};


	auto partial_m_K = [&](int a, int b, int dim){
        if(dim==0) return (globalGrid[i+1][j][k].K[a][b] - globalGrid[i-1][j][k].K[a][b])/(2.*DX);
        if(dim==1) return (globalGrid[i][j+1][k].K[a][b] - globalGrid[i][j-1][k].K[a][b])/(2.*DY);
        if(dim==2) return (globalGrid[i][j][k+1].K[a][b] - globalGrid[i][j][k-1].K[a][b])/(2.*DZ);
        return 0.0;
    };

    auto partial_m_gamma = [&](int a, int b, int dim){
        if(dim==0) return (globalGrid[i+1][j][k].gamma[a][b] - globalGrid[i-1][j][k].gamma[a][b])/(2.*DX);
        if(dim==1) return (globalGrid[i][j+1][k].gamma[a][b] - globalGrid[i][j-1][k].gamma[a][b])/(2.*DY);
        if(dim==2) return (globalGrid[i][j][k+1].gamma[a][b] - globalGrid[i][j][k-1].gamma[a][b])/(2.*DZ);
        return 0.0;
    };

    auto partial_j_beta_comp = [&](int dimSpace, int compBeta){
        if(dimSpace==0) return partialX_betacomp(grid_obj, i, j, k, compBeta);
        if(dimSpace==1) return partialY_betacomp(grid_obj, i, j, k, compBeta);
        if(dimSpace==2) return partialZ_betacomp(grid_obj, i, j, k, compBeta);
        return 0.0;
    };

    double dtGamma[3][3];
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            double rhs = -2.0 * alpha * KLocal[a][b];

            double advection = 0.0;
            for(int m=0; m<3; m++){
                advection += beta[m] * partial_m_gamma(a,b,m);
            }


			double shiftTerm = 0.0;
			for (int m = 0; m < 3; m++) {
				shiftTerm += gammaLocal[a][m] * partial_j_beta_comp(m, b);
				shiftTerm += gammaLocal[b][m] * partial_j_beta_comp(m, a);
			}


            rhs += advection + shiftTerm;

            dtGamma[a][b] = rhs;
        }
    }

    double dtK[3][3];
	for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            double dDaDb_alpha = - D_iD_j_alpha(a,b);
            double termB = 0.0; 
            {
                double KContra[3][3];
				#pragma omp simd collapse(2)
				for(int m=0; m<3; m++){
                    for(int n=0; n<3; n++){
                        double sum_ = 0.0;
                        for(int r=0; r<3; r++){
                            for(int s=0; s<3; s++){
                                sum_ += gammaInv[m][r]*gammaInv[n][s]*KLocal[r][s];
                            }
                        }
                        KContra[m][n] = sum_;
                    }
                }
                double KtraceLocal = 0.0;
                for(int x=0; x<3; x++){
                    for(int y=0; y<3; y++){
                        KtraceLocal += gammaInv[x][y]*KLocal[x][y];
                    }
                }
                double RicciTerm = Ricci[a][b];
                double KKab      = KtraceLocal * KLocal[a][b];
                
                double KaM_KM_b = 0.0;
                for(int m=0; m<3; m++){
                    KaM_KM_b += KLocal[a][m]*KContra[m][b];
                }
                termB = alpha*(RicciTerm + KKab - 2.0*KaM_KM_b);
            }

            double lieBetaK = 0.0;
            {
                double advK = 0.0;
                for(int m=0; m<3; m++){
                    advK += beta[m]*partial_m_K(a,b,m);
                }
                double shiftK = 0.0;
                for(int m=0; m<3; m++){
                    shiftK += KLocal[m][b] * partial_j_beta_comp(a,m);
                    shiftK += KLocal[a][m] * partial_j_beta_comp(b,m);
                }
                lieBetaK = advK + shiftK;
            }

            dtK[a][b] = dDaDb_alpha + termB + lieBetaK;
        }
    }

    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            cell.dgt[a][b] = dtGamma[a][b];
            cell.dKt[a][b] = dtK[a][b];
        }
    }

	double norm_dtK = 0.0;
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			norm_dtK += dtK[a][b] * dtK[a][b];
		}
	}
	norm_dtK = sqrt(norm_dtK / 9.0); 


	/* printf(" norm(dtK) = %e\n", norm_dtK); */
}



