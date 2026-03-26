#pragma once

#include "SimdConfig.h"
#include "core/Constants.h"

#include <cstddef>

typedef struct {
  float x, y, z;
  float lambda;
} GeodesicPoint;

using VEC_TYPE = curvatureengine::simd::Vec4d;

using ChristoffelTensor = double[NDIM][NDIM][NDIM];
using ChristoffelEvalFn = void (*)(const double coords[NDIM],
                                   ChristoffelTensor gamma, void *ctx);

using ChristoffelTensorVec = VEC_TYPE[NDIM][NDIM][NDIM];
using ChristoffelEvalFnVec = void (*)(const VEC_TYPE coords[NDIM],
                                      ChristoffelTensorVec &gamma, void *ctx);

void geodesic_AVX(VEC_TYPE x[NDIM], VEC_TYPE v[NDIM], float lambda_max,
                  ChristoffelEvalFnVec evaluator, void *ctx,
                  VEC_TYPE step_size);
VEC_TYPE geodesic_raytrace_AVX(VEC_TYPE x[NDIM], VEC_TYPE v[NDIM],
                               float lambda_max, ChristoffelEvalFnVec evaluator,
                               void *ctx, VEC_TYPE step_size);
bool initialize_geodesic_point_storage(std::size_t max_points);
int get_stored_geodesic_point_count();
void release_geodesic_point_storage();
void store_geodesic_point(float x[4], float lambda);
void store_geodesic_point_AVX(VEC_TYPE x[NDIM], float lambda);
void write_vtk_file(const char *filename);
