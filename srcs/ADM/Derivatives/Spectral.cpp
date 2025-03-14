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
