#include <Geodesics.h>
#include <Derivatives.h>
void Grid::compute_spectral_derivatives_for_gamma() {
    Derivatives deriv;
	for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            std::vector<double> gamma_ab(NX*NY*NZ), dgammaX_ab(NX*NY*NZ), dgammaY_ab(NX*NY*NZ), dgammaZ_ab(NX*NY*NZ);
	
            for (int ix = 0; ix < NX; ix++) {
                for (int iy = 0; iy < NY; iy++) {
                    for (int iz = 0; iz < NZ; iz++) {
                        int index = (ix*(NY*NZ)) + (iy*NZ) + iz;
                        gamma_ab[index] = this->getCell(ix, iy, iz).gamma[a][b];
                    }
                }
            }

            spectral_derivative_3D(gamma_ab, dgammaX_ab, NX, NY, NZ, DX, DY, DZ, 0); 
            spectral_derivative_3D(gamma_ab, dgammaY_ab, NX, NY, NZ, DX, DY, DZ, 1); 
            spectral_derivative_3D(gamma_ab, dgammaZ_ab, NX, NY, NZ, DX, DY, DZ, 2);

            this->dgammaX[a][b] = dgammaX_ab;
            this->dgammaY[a][b] = dgammaY_ab;
            this->dgammaZ[a][b] = dgammaZ_ab;
        }
    }
}

