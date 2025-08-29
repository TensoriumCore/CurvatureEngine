#include <Geodesics.h>


void BSSNevolve::compute_dt_chi(Grid &grid_obj, int i, int j, int k, float &dt_chi) {
    const auto &cell = grid_obj.getCell(i, j, k);

    const float chi   = cell.chi;
    const float alpha = cell.gauge.alpha;

    float Ktrace = 0.0f;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        Ktrace += cell.geom.tildgamma_inv[a][b] * cell.curv.K[a][b];
    Ktrace *= chi;

    float div_beta = 0.0f;
    for (int a = 0; a < 3; ++a) {
        div_beta += partial_m(grid_obj, i, j, k, a, [&](const Grid::Cell2D &c){
            return c.gauge.beta[a];
        });
    }

    float beta_grad_chi = 0.0f;
    for (int a = 0; a < 3; ++a) {
        float dchi = partial_m(grid_obj, i, j, k, a, [&](const Grid::Cell2D &c){
            return c.chi;
        });
        beta_grad_chi += cell.gauge.beta[a] * dchi;
    }

    float T1 = (2.0f/3.0f) * chi * alpha * Ktrace; 
    float T2 = -(2.0f/3.0f) * chi * div_beta; 
    float T3 = beta_grad_chi;
    dt_chi = T1 + T2 + T3;

    // S := (∂t χ − β·∇χ)/χ  = (2/3)(α K − ∂_i β^i)
    float S = (dt_chi - T3) / (chi + 1e-30f);
    float K_meas = ( (1.5f)*S + div_beta ) / (alpha + 1e-30f); 

}



