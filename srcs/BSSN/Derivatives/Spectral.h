#pragma once

#include <Geodesics.h>

class ChebyshevSpectral {
public:
    static std::vector<float> chebyshev_nodes(int N, float a, float b) {
        std::vector<float> nodes(N);
        for (int k = 0; k < N; ++k) {
            nodes[k] = 0.5 * ((b - a) * cos(M_PI * k / (N - 1)) + (b + a));
        }
        return nodes;
    }

    static std::vector<std::vector<float>> chebyshev_diff_matrix(const std::vector<float>& nodes, float a, float b) {
        int N = nodes.size();
        std::vector<std::vector<float>> D(N, std::vector<float>(N, 0.0));
        std::vector<float> c(N, 1.0);
        c[0] = 2.0;
        c[N - 1] = 2.0;
        float factor = 2.0 / (b - a);
        
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i != j) {
                    D[i][j] = factor * (c[i] / c[j]) * pow(-1, i + j) / (nodes[i] - nodes[j]);
                }
            }
        }
        for (int i = 0; i < N; ++i) {
            float sum = 0.0;
            for (int j = 0; j < N; ++j) {
                if (i != j)
                    sum += D[i][j];
            }
            D[i][i] = -sum;
        }
        return D;
    }

    static std::vector<float> spectral_derivative(const std::vector<float>& f, const std::vector<std::vector<float>>& D) {
        int N = f.size();
        std::vector<float> df(N, 0.0);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                df[i] += D[i][j] * f[j];
            }
        }
        return df;
    }
};
