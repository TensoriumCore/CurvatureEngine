#pragma once

#include "core/Constants.h"

#include <array>

using Matrix2x2 = std::array<std::array<float, 2>, 2>;
using Matrix3x3 = std::array<std::array<float, DIM3>, DIM3>;
using Matrix4x4 = std::array<std::array<float, NDIM>, NDIM>;
using MatrixNDIM = std::array<std::array<float, NDIM>, NDIM>;
using Vector3 = std::array<float, DIM3>;
using Vector4 = std::array<float, NDIM>;
using Tensor3D =
    std::array<std::array<std::array<float, DIM3>, DIM3>, DIM3>;
using Tensor4D = std::array<
    std::array<std::array<std::array<float, DIM3>, DIM3>, DIM3>, DIM3>;
using Christoffel3D =
    std::array<std::array<std::array<float, DIM3>, DIM3>, DIM3>;
using Riemann3D = std::array<
    std::array<std::array<std::array<float, DIM3>, DIM3>, DIM3>, DIM3>;
