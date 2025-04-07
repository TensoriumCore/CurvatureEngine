#include <Geodesics.h>
#include <cassert>


void Grid::compute_time_derivatives(Grid &grid_obj, int i, int j, int k)
{
    Cell2D &cell = globalGrid[i][j][k];
    BSSNevolve bssn;
    float alpha = cell.gauge.alpha;
    float beta[3] = { cell.gauge.beta[0], cell.gauge.beta[1], cell.gauge.beta[2] };

      GridTensor gridTensor;
    float Gamma[3][3][3];
    gridTensor.compute_christoffel_3D(grid_obj, i, j, k, Gamma);

    float Ricci[3][3];
    gridTensor.compute_ricci_BSSN(grid_obj, i, j, k, Ricci);

    float partialBeta[3][3];
    for (int dim = 0; dim < 3; ++dim) {
        for (int comp = 0; comp < 3; ++comp) {
            partialBeta[dim][comp] = partial_m(grid_obj, i, j, k, dim,
                [&](const Grid::Cell2D &c) { return c.gauge.beta[comp]; }
            );
			cell.gauge.dt_beta[comp] = partialBeta[dim][comp];
        }
    }


	float partialAlpha[3];
	for (int dim = 0; dim < 3; ++dim) {
		partialAlpha[dim] = partial_m(grid_obj, i, j, k, dim,
				[&](const Grid::Cell2D &c) { return c.gauge.alpha; }
				);
	}

	float f_alpha = 1.0; 
	cell.gauge.dt_alpha = -alpha * alpha * f_alpha * cell.curv.K_trace;
	for (int m = 0; m < 3; ++m) {
		cell.gauge.dt_alpha += beta[m] * partialAlpha[m];
	}


    float d2Alpha[3][3];
    for (int m2 = 0; m2 < 3; ++m2) {
        for (int n2 = 0; n2 < 3; ++n2) {
			d2Alpha[m2][n2] = second_partial_alpha(grid_obj, i, j, k, m2, n2);
		}
    }
    float partialKtrace[3];
    for (int dim = 0; dim < 3; ++dim) {
        partialKtrace[dim] = partial_m(grid_obj, i, j, k, dim,
            [&](const Grid::Cell2D &c) { return c.curv.K_trace; }
        );
    }

    float partialTildeGamma[3][3][3];
    float partialAtilde[3][3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int dim = 0; dim < 3; dim++) {
                partialTildeGamma[dim][a][b] = partial_m(grid_obj, i, j, k, dim,
                    [&](const Grid::Cell2D &c) { return c.geom.tilde_gamma[a][b]; }
                );
                partialAtilde[dim][a][b] = partial_m(grid_obj, i, j, k, dim,
                    [&](const Grid::Cell2D &c) { return c.atilde.Atilde[a][b]; }
                );
            }
        }
    }

    float dt_chi = 0.0;
    bssn.compute_dt_chi(grid_obj, i, j, k, dt_chi);
    cell.dt_chi = dt_chi;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            float adv = 0.0;
            for (int m = 0; m < 3; ++m) {
                adv += beta[m] * partialTildeGamma[m][a][b];
            }

            float shift = 0.0;
            for (int m = 0; m < 3; ++m) {
                shift += cell.geom.tilde_gamma[a][m] * partialBeta[m][b];
                shift += cell.geom.tilde_gamma[b][m] * partialBeta[m][a];
            }

			float div_beta = 0.0;
			for (int m = 0; m < 3; ++m) {
				div_beta += partialBeta[m][m]; 
				for (int l = 0; l < 3; ++l) {
					div_beta += Gamma[m][m][l] * beta[l];  
				}
			}

			cell.geom.dt_tilde_gamma[a][b] = -2.0 * alpha * cell.atilde.Atilde[a][b] 
				+ adv 
				+ shift 
				- (2.0/3.0) * cell.geom.tilde_gamma[a][b] * div_beta;
			cell.dgt[a][b] = cell.geom.dt_tilde_gamma[a][b];
        }
    }
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {

            float D2_alpha = d2Alpha[a][b];
            float sumG = 0.0;
            for (int m = 0; m < 3; ++m) {
                sumG += Gamma[m][a][b] * partialAlpha[m];
            }
            D2_alpha -= sumG;

            float trace_D2_alpha = 0.0;
#pragma omp simd reduction(+:trace_D2_alpha)
            for (int mm = 0; mm < 3; ++mm) {
                for (int nn = 0; nn < 3; ++nn) {
                    float part = d2Alpha[mm][nn];
                    float gpart = 0.0;
                    for (int ll = 0; ll < 3; ++ll) {
                        gpart += Gamma[ll][mm][nn] * partialAlpha[ll];
                    }
                    trace_D2_alpha += cell.geom.tildgamma_inv[mm][nn] * (part - gpart);
                }
            }

            float Ricci_TF = Ricci[a][b];
            float R_scalar = 0.0;
            for (int mm = 0; mm < 3; ++mm) {
                for (int nn = 0; nn < 3; ++nn) {
                    R_scalar += cell.geom.tildgamma_inv[mm][nn] * Ricci[mm][nn];
                }
            }
            Ricci_TF -= (1.0/3.0) * cell.geom.tilde_gamma[a][b] * R_scalar;

            float A_A = 0.0;
            for (int k1 = 0; k1 < 3; ++k1) {
                for (int l1 = 0; l1 < 3; ++l1) {
                    A_A += cell.atilde.Atilde[a][k1]
                         * cell.geom.tildgamma_inv[k1][l1]
                         * cell.atilde.Atilde[l1][b];
                }
            }

            float adv = 0.0;
            for (int m = 0; m < 3; ++m) {
                adv += beta[m] * partialAtilde[m][a][b];
            }

            float shift_term = 0.0;
            for (int m = 0; m < 3; ++m) {
                shift_term += cell.atilde.Atilde[m][b] * partialBeta[m][a];
                shift_term += cell.atilde.Atilde[a][m] * partialBeta[m][b];
            }

            cell.atilde.dt_Atilde[a][b] =
                  cell.chi * (
                      -D2_alpha
                      + (1.0/3.0) * cell.geom.tilde_gamma[a][b] * trace_D2_alpha
                      + alpha * Ricci_TF
                  )
                + alpha * (
                      cell.curv.K_trace * cell.atilde.Atilde[a][b]
                      - 2.0 * A_A
                  )
                + adv
                + shift_term;
        }
    }
    float laplacian_alpha = 0.0;
    for (int mm = 0; mm < 3; ++mm) {
        for (int nn = 0; nn < 3; ++nn) {
            float d2a = d2Alpha[mm][nn];
            float gamma_conn = 0.0;
#pragma omp simd
            for (int ll = 0; ll < 3; ++ll) {
                gamma_conn += Gamma[ll][mm][nn] * partialAlpha[ll];
            }
            laplacian_alpha += cell.geom.tildgamma_inv[mm][nn] * (d2a - gamma_conn);
        }
    }

	float Atilde_raised[3][3] = {0.0};
#pragma omp parallel for collapse(2) 
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l) {
					Atilde_raised[i][j] += cell.geom.tildgamma_inv[i][k]
						* cell.geom.tildgamma_inv[j][l]
						* cell.atilde.Atilde[k][l];
				}
			}
		}
	}

	float Atilde_squared = 0.0;
#pragma omp simd reduction(+:Atilde_squared)  
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			Atilde_squared += cell.atilde.Atilde[i][j] * Atilde_raised[i][j];
		}
	}
	float total_R_scalar = 0.0;
#pragma omp simd reduction(+:total_R_scalar)
    for (int a1 = 0; a1 < 3; ++a1) {
        for (int b1 = 0; b1 < 3; ++b1) {
            total_R_scalar += cell.geom.tildgamma_inv[a1][b1] * Ricci[a1][b1];
        }
    }

    float adv_K = 0.0;
    for (int m = 0; m < 3; ++m) {
        adv_K += beta[m] * partialKtrace[m];
    }

    cell.curv.dt_K_trace =
        -laplacian_alpha
        + alpha * (
              Atilde_squared
            + (1.0/3.0)*cell.curv.K_trace*cell.curv.K_trace
            + total_R_scalar
          )
        + adv_K;
}
