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

/* void spectral_derivative_3D(const std::vector<float> &f_in, */
/*                             std::vector<float> &f_out, */
/*                             int N_x, int N_y, int N_z, */
/*                             float dx, float dy, float dz, */
/*                             int dim) */
/* { */
/*     int N = N_x * N_y * N_z; */
/*     fftw_complex *in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); */
/*     fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); */
/*  */
/*     for (int i = 0; i < N; i++) { */
/*         in[i][0] = f_in[i]; // réel */
/*         in[i][1] = 0.0;     // imag */
/*     } */
/*  */
/*     fftw_plan plan_forward  = fftw_plan_dft_3d(N_x, N_y, N_z, in, out, FFTW_FORWARD,  FFTW_ESTIMATE); */
/*     fftw_plan plan_backward = fftw_plan_dft_3d(N_x, N_y, N_z, out, in, FFTW_BACKWARD, FFTW_ESTIMATE); */
/*  */
/*     // FFT directe 3D */
/*     fftw_execute(plan_forward); */
/*  */
/*     for (int ix = 0; ix < N_x; ix++) { */
/*         int kx = (ix <= N_x/2) ? ix : ix - N_x; */
/*         float fx = 2.0 * M_PI * kx / (N_x * dx); */
/*  */
/*         for (int iy = 0; iy < N_y; iy++) { */
/*             int ky = (iy <= N_y/2) ? iy : iy - N_y; */
/*             float fy = 2.0 * M_PI * ky / (N_y * dy); */
/*  */
/*             for (int iz = 0; iz < N_z; iz++) { */
/*                 int kz = (iz <= N_z/2) ? iz : iz - N_z; */
/*                 float fz = 2.0 * M_PI * kz / (N_z * dz); */
/*  */
/*                 int index = (ix*(N_y*N_z)) + (iy*N_z) + iz; */
/*                 float real = out[index][0]; */
/*                 float imag = out[index][1]; */
/*  */
/*                 float factor = 0.0; */
/*                 if      (dim == 0) factor = fx;   */
/*                 else if (dim == 1) factor = fy;   */
/*                 else if (dim == 2) factor = fz;   */
/*  */
/*                 out[index][0] = -factor * imag; */
/*                 out[index][1] =  factor * real;  */
/*             } */
/*         } */
/*     } */
/*  */
/*     fftw_execute(plan_backward); */
/*  */
/*     for (int i = 0; i < N; i++) { */
/*         f_out[i] = in[i][0] / float(N); */
/*     } */
/*  */
/*     fftw_destroy_plan(plan_forward); */
/*     fftw_destroy_plan(plan_backward); */
/*     fftw_free(in); */
/*     fftw_free(out); */
/* } */
/*  */
/* std::vector<float> computeBaryWeights(const std::vector<float>& nodes) { */
/*     int N = nodes.size(); */
/*     std::vector<float> weights(N, 1.0); */
/*     for (int j = 0; j < N; j++) { */
/*         weights[j] = 1.0; */
/*         for (int k = 0; k < N; k++) { */
/*             if (k != j) { */
/*                 weights[j] /= (nodes[j] - nodes[k]); */
/*             } */
/*         } */
/*     } */
/*     return weights; */
/* } */
/*  */
/*
 * This function compute the barycentric interpolation of a function
 * @param x the point to interpolate
 * @param nodes the nodes of the function
 * @param fVals the values of the function at the nodes
 * @param weights the weights of the barycentric interpolation
 * @return the interpolated value
 * */

/* float barycentricInterpolate(float x,  */
/*                                 const std::vector<float>& nodes,  */
/*                                 const std::vector<float>& fVals,  */
/*                                 const std::vector<float>& weights) { */
/*     const float tol = 1e-12; */
/*     int N = nodes.size(); */
/*     float numerator = 0.0; */
/*     float denominator = 0.0; */
/*     for (int j = 0; j < N; j++) { */
/*         if (std::fabs(x - nodes[j]) < tol) { */
/*             return fVals[j]; */
/*         } */
/*         float temp = weights[j] / (x - nodes[j]); */
/*         numerator += temp * fVals[j]; */
/*         denominator += temp; */
/*     } */
/*     return numerator / denominator; */
/* } */
/*  */
/* std::vector<float> interpolateUniformToNodes_barycentric(const std::vector<float>& xUniform, */
/*                                                           const std::vector<float>& fUniform, */
/*                                                           const std::vector<float>& xTarget) { */
/*     std::vector<float> weights = computeBaryWeights(xUniform); */
/*     std::vector<float> fTarget; */
/*     fTarget.reserve(xTarget.size()); */
/*     for (float x : xTarget) { */
/*         fTarget.push_back(barycentricInterpolate(x, xUniform, fUniform, weights)); */
/*     } */
/*     return fTarget; */
/* } */
/*  */
/* float GridTensor::partialX_gammaSpec(Grid &grid_obj, int i, int j, int k, int a, int b) { */
/*     std::vector<float> gamma_vals; */
/*     for (int idx = 0; idx < NX; ++idx) { */
/*         gamma_vals.push_back(grid_obj.getCell(idx, j, k).gamma[a][b]); */
/*     } */
/*      */
/*     std::vector<float> xUniform; */
/*     float x_min = -128.0; */
/*     float x_max = 128.0; */
/*     for (int idx = 0; idx < NX; ++idx) { */
/*         float x = x_min + idx * (x_max - x_min) / (NX - 1); */
/*         xUniform.push_back(x); */
/*     } */
/*      */
/*     std::vector<float> chebNodes = ChebyshevSpectral::chebyshev_nodes(NX, x_min, x_max); */
/*     std::vector<float> gamma_vals_nodes = interpolateUniformToNodes_barycentric(xUniform, gamma_vals, chebNodes); */
/*     auto D = ChebyshevSpectral::chebyshev_diff_matrix(chebNodes, x_min, x_max); */
/*     std::vector<float> d_gamma_nodes = ChebyshevSpectral::spectral_derivative(gamma_vals_nodes, D); */
/*     std::vector<float> d_gamma = interpolateUniformToNodes_barycentric(chebNodes, d_gamma_nodes, xUniform); */
/*     return d_gamma[i]; */
/* } */
/*  */
