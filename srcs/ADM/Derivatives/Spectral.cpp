#include "Spectral.h"

/*
 * This function compute Partial derivative using spectral methods
 * its based on the Chebyshev spectral methods and the barycentric interpolation
 *
 * This is not currently working because the Grid is not currently ready to use theses differentiation methods
 * @param i , j , k the index of the cell
 * @param a , b the index of the tensor
 * @return the partial derivative of the tensor
 *
 * THIS IS JUST AN EXAMPLE OF HOW TO USE THE SPECTRAL METHODS
 *
 * */

void spectral_derivative_3D(const std::vector<double> &f_in,
                            std::vector<double> &f_out,
                            int N_x, int N_y, int N_z,
                            double dx, double dy, double dz,
                            int dim)
{
    int N = N_x * N_y * N_z;
    fftw_complex *in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        in[i][0] = f_in[i]; // réel
        in[i][1] = 0.0;     // imag
    }

    fftw_plan plan_forward  = fftw_plan_dft_3d(N_x, N_y, N_z, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
    fftw_plan plan_backward = fftw_plan_dft_3d(N_x, N_y, N_z, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    // FFT directe 3D
    fftw_execute(plan_forward);

    for (int ix = 0; ix < N_x; ix++) {
        int kx = (ix <= N_x/2) ? ix : ix - N_x;
        double fx = 2.0 * M_PI * kx / (N_x * dx);

        for (int iy = 0; iy < N_y; iy++) {
            int ky = (iy <= N_y/2) ? iy : iy - N_y;
            double fy = 2.0 * M_PI * ky / (N_y * dy);

            for (int iz = 0; iz < N_z; iz++) {
                int kz = (iz <= N_z/2) ? iz : iz - N_z;
                double fz = 2.0 * M_PI * kz / (N_z * dz);

                int index = (ix*(N_y*N_z)) + (iy*N_z) + iz;
                double real = out[index][0];
                double imag = out[index][1];

                double factor = 0.0;
                if      (dim == 0) factor = fx;  
                else if (dim == 1) factor = fy;  
                else if (dim == 2) factor = fz;  

                out[index][0] = -factor * imag;
                out[index][1] =  factor * real; 
            }
        }
    }

    fftw_execute(plan_backward);

    for (int i = 0; i < N; i++) {
        f_out[i] = in[i][0] / double(N);
    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);
}

std::vector<double> computeBaryWeights(const std::vector<double>& nodes) {
    int N = nodes.size();
    std::vector<double> weights(N, 1.0);
    for (int j = 0; j < N; j++) {
        weights[j] = 1.0;
        for (int k = 0; k < N; k++) {
            if (k != j) {
                weights[j] /= (nodes[j] - nodes[k]);
            }
        }
    }
    return weights;
}

/*
 * This function compute the barycentric interpolation of a function
 * @param x the point to interpolate
 * @param nodes the nodes of the function
 * @param fVals the values of the function at the nodes
 * @param weights the weights of the barycentric interpolation
 * @return the interpolated value
 * */

double barycentricInterpolate(double x, 
                                const std::vector<double>& nodes, 
                                const std::vector<double>& fVals, 
                                const std::vector<double>& weights) {
    const double tol = 1e-12;
    int N = nodes.size();
    double numerator = 0.0;
    double denominator = 0.0;
    for (int j = 0; j < N; j++) {
        if (std::fabs(x - nodes[j]) < tol) {
            return fVals[j];
        }
        double temp = weights[j] / (x - nodes[j]);
        numerator += temp * fVals[j];
        denominator += temp;
    }
    return numerator / denominator;
}

std::vector<double> interpolateUniformToNodes_barycentric(const std::vector<double>& xUniform,
                                                          const std::vector<double>& fUniform,
                                                          const std::vector<double>& xTarget) {
    std::vector<double> weights = computeBaryWeights(xUniform);
    std::vector<double> fTarget;
    fTarget.reserve(xTarget.size());
    for (double x : xTarget) {
        fTarget.push_back(barycentricInterpolate(x, xUniform, fUniform, weights));
    }
    return fTarget;
}

double GridTensor::partialX_gammaSpec(Grid &grid_obj, int i, int j, int k, int a, int b) {
    std::vector<double> gamma_vals;
    for (int idx = 0; idx < NX; ++idx) {
        gamma_vals.push_back(grid_obj.getCell(idx, j, k).gamma[a][b]);
    }
    
    std::vector<double> xUniform;
    double x_min = -128.0;
    double x_max = 128.0;
    for (int idx = 0; idx < NX; ++idx) {
        double x = x_min + idx * (x_max - x_min) / (NX - 1);
        xUniform.push_back(x);
    }
    
    std::vector<double> chebNodes = ChebyshevSpectral::chebyshev_nodes(NX, x_min, x_max);
    std::vector<double> gamma_vals_nodes = interpolateUniformToNodes_barycentric(xUniform, gamma_vals, chebNodes);
    auto D = ChebyshevSpectral::chebyshev_diff_matrix(chebNodes, x_min, x_max);
    std::vector<double> d_gamma_nodes = ChebyshevSpectral::spectral_derivative(gamma_vals_nodes, D);
    std::vector<double> d_gamma = interpolateUniformToNodes_barycentric(chebNodes, d_gamma_nodes, xUniform);
    return d_gamma[i];
}

