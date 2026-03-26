#include "Connexion.h"
#include "Metric.h"
#include "Tensor.h"
#include "app/RuntimeState.h"
#include "core/Constants.h"

#include <cmath>
#include <iostream>

namespace {

constexpr float kRiemannDelta = 5.0e-2f;
constexpr float kPi = 3.14159265358979323846f;

struct CurvatureSummary {
    Tensor::Riemann4D riemann{};
    Tensor::MatrixNDIM ricci{};
    double ricci_scalar = 0.0;
    double kretschmann = 0.0;
};

CurvatureSummary compute_curvature(const char* metric_name, float spin, float radius) {
    a = spin;

    Metric metric;
    Connexion connexion;
    Tensor tensor;
    CurvatureSummary summary;
    Connexion::Christoffel3D gamma{};
    Connexion::Christoffel4D gamma_plus_h{};
    Connexion::Christoffel4D gamma_minus_h{};
    std::array<float, NDIM> x = {0.0f, radius, kPi / 2.0f, 0.0f};

    metric.calculate_metric_by_name(metric_name, x, metric.gcov, metric.gcon);
    connexion.calculate_christoffel(x, kRiemannDelta, gamma, metric.gcov, metric.gcon, metric_name);

    for (int d = 0; d < NDIM; ++d) {
        tensor.calculate_Gamma_at_offset(x, d, kRiemannDelta, kRiemannDelta,
                                         metric.gcov, metric.gcon, gamma_plus_h[d],
                                         metric_name);
        tensor.calculate_Gamma_at_offset(x, d, -kRiemannDelta, kRiemannDelta,
                                         metric.gcov, metric.gcon, gamma_minus_h[d],
                                         metric_name);
    }

    tensor.calculate_riemann(gamma, gamma_plus_h, gamma_minus_h, summary.riemann,
                             kRiemannDelta);
    tensor.contract_riemann(summary.riemann, summary.ricci);
    summary.ricci_scalar = tensor.calculate_ricci_scalar(summary.ricci, metric.gcon);
    summary.kretschmann =
        tensor.calculate_kretschmann_scalar(summary.riemann, metric.gcov, metric.gcon);

    return summary;
}

bool expect_near(const char* label, double actual, double expected, double tolerance) {
    const double error = std::fabs(actual - expected);
    if (error <= tolerance) {
        return true;
    }

    std::cerr << label << " mismatch: got " << actual << ", expected " << expected
              << " +/- " << tolerance << '\n';
    return false;
}

bool expect_positive(const char* label, double actual) {
    if (actual > 0.0) {
        return true;
    }

    std::cerr << label << " expected positive value, got " << actual << '\n';
    return false;
}

}  // namespace

int main() {
    bool ok = true;

    {
        const auto summary = compute_curvature("schwarzschild", 0.0f, 6.0f);
        const double expected_kretschmann = 48.0 / std::pow(6.0, 6.0);
        ok &= expect_near("Schwarzschild Ricci scalar", summary.ricci_scalar, 0.0, 1.0e-3);
        ok &= expect_near("Schwarzschild Kretschmann", summary.kretschmann,
                          expected_kretschmann, 5.0e-5);
        ok &= expect_near("Schwarzschild R^0_{101}",
                          static_cast<double>(summary.riemann[0][1][0][1]),
                          1.0 / 72.0, 5.0e-3);
        ok &= expect_near("Schwarzschild R^1_{212}",
                          static_cast<double>(summary.riemann[1][2][1][2]),
                          -1.0 / 6.0, 5.0e-3);
    }

    {
        const auto summary = compute_curvature("kerr-newman", 0.935f, 6.0f);
        ok &= expect_near("Kerr-Newman Ricci scalar", summary.ricci_scalar, 0.0,
                          1.0e-3);
        ok &= expect_positive("Kerr-Newman Kretschmann", summary.kretschmann);
    }

    {
        const auto summary = compute_curvature("kds", 0.35f, 6.0f);
        ok &= expect_near("Kerr-de Sitter Ricci scalar", summary.ricci_scalar,
                          4.0 * KERR_DE_SITTER_LAMBDA, 2.0e-4);
        ok &= expect_positive("Kerr-de Sitter Kretschmann", summary.kretschmann);
    }

    return ok ? 0 : 1;
}
