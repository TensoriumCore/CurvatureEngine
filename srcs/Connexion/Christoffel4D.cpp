#include <Geodesics.h>


void Connexion::calculate_christoffel(const VectorNDIM& X, float h,
                                      Christoffel3D& gamma,
                                      std::array<std::array<float, NDIM>, NDIM>& g,
                                      std::array<std::array<float, NDIM>, NDIM>& g_inv, 
                                      const char* metric) {
    Christoffel3D tmp{};
    VectorNDIM Xh = X, Xl = X;
    MatrixNDIM gh{}, gl{};
    Metric metric_obj;

    gamma.fill({}); 

    for (int mu = 0; mu < NDIM; mu++) {
        Xh = X;
        Xl = X;
        Xh[mu] += DELTA;
        Xl[mu] -= DELTA;

        if (strcmp(metric, "schwarzschild") == 0 || 
            strcmp(metric, "kerr") == 0 ||
            strcmp(metric, "kerr-newman") == 0 ||
            strcmp(metric, "ds") == 0) {
            
			metric_obj.calculate_metric(Xh, gh, g_inv); 
            metric_obj.calculate_metric(Xl, gl, g_inv);
        }

        for (int lam = 0; lam < NDIM; lam++) {
            for (int nu = 0; nu < NDIM; nu++) {
                gamma[lam][nu][mu] = (gh[lam][nu] - gl[lam][nu]) / (2 * DELTA);
            }
        }
    }

    for (int lam = 0; lam < NDIM; lam++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int mu = 0; mu < NDIM; mu++) {
                tmp[lam][nu][mu] = 0.5 * (gamma[nu][lam][mu] + 
                                          gamma[mu][lam][nu] - 
                                          gamma[mu][nu][lam]);
            }
        }
    }

    for (int lam = 0; lam < NDIM; lam++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int mu = 0; mu < NDIM; mu++) {
                gamma[lam][nu][mu] = 0.0;
                for (int kap = 0; kap < NDIM; kap++) {
                    gamma[lam][nu][mu] += g_inv[lam][kap] * tmp[kap][nu][mu];
                }
            }
        }
    }    

    std::cout << "Christoffel symbols calculated\n";
    print_christoffel_matrix(gamma); 
    check_symmetry_christoffel(gamma);
}


