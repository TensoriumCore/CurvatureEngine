
#pragma once

#include "SimdConfig.h"
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <vector>
// #include <fftw3.h>
#define C 1.0
#define G 1.0
#define M 1.0
#define BLOCK_SIZE 1024
#define BUFFER_SIZE 1024
#define SMALL 1.e-40
#define NDIM 4
#define DT 0.0000005
#define max_dt 19999.0
#define ALIGNMENT 32
#if CURVATUREENGINE_TARGET_AVX2
#define ARCH "AVX2"
#elif CURVATUREENGINE_TARGET_NEON
#define ARCH "NEON"
#else
#define ARCH "SCALAR"
#endif
#define TOLERANCE 1e-10
#define DELTA 1e-6
#define NDIM3 3
#define DELTA3 1e-4

using Matrix2x2 = std::array<std::array<float, 2>, 2>;
using Matrix3x3 = std::array<std::array<float, 3>, 3>;
using Matrix4x4 = std::array<std::array<float, 4>, 4>;
using MatrixNDIM = std::array<std::array<float, NDIM>, NDIM>;

#include <Connexion.h>
#include <Derivatives.h>
#include <Grid.h>
#include <GridTensor.h>
#include <Log.h>
#include <Metric.h>
#include <Tensor.h>
#include <matrix.h>

typedef struct {
  float x, y, z;
  float lambda;
} GeodesicPoint;

using VEC_TYPE = curvatureengine::simd::Vec4d;

using ChristoffelTensor = double[NDIM][NDIM][NDIM];
using ChristoffelEvalFn = void (*)(const double coords[NDIM],
                                   ChristoffelTensor gamma, void *ctx);

void write_vtk_file(const char *filename);
void store_geodesic_point(float x[4], float lambda);
void geodesic_AVX(VEC_TYPE x[4], VEC_TYPE v[4], float lambda_max,
                  ChristoffelEvalFn evaluator, void *ctx, VEC_TYPE step_size);
void store_geodesic_point_AVX(VEC_TYPE x[4], float lambda);

float calculate_impact_parameter(float p_t, float p_phi, float g_tt,
                                 float g_tphi, float g_phiphi);
float calculate_emission_angle(float p_r, float p_phi, float g_rr,
                               float g_phiphi);
float b_critique_kerr(float a, int sense);
int compute_photon_properties(float g[4][4], float p[4]);

/* Problem specific functions */

int Riemann_tensor(const char *metric);
int Geodesics_prob();
int light_geodesics_prob();
int Metric_prob();
int grid_setup();
void generate_blackhole_image();
void generate_blackhole_shadow();
void evolveADM(Grid::Cell2D &cell, int i, int j, float dt,
               const std::vector<std::vector<float>> &alpha_grid, float r_min,
               float theta_min, float dr, float dtheta, std::ofstream &file,
               int step);
void calc_gamma_ij_2D(int i, int j, float r_min, float dr, float theta_min,
                      float dtheta, Metric &metric_obj, Matrix3x3 &gamma3,
                      Matrix3x3 &gamma3_inv);
void calc_gamma_ij(const Vector3 &X3D, Matrix3x3 &gamma3,
                   Matrix3x3 &gamma3_inv);
