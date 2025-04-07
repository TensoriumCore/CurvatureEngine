#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <time.h>
#include <chrono>
#include <sys/time.h>
#include <iostream>
#include <array>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>
#include <fftw3.h>
#define C 1.0
#define G 1.0 
#define M 1.0
#define BLOCK_SIZE 1024
#define BUFFER_SIZE 1024
#define SMALL 1.e-40
#define NDIM 4
#define DT 0.0000005
#define max_dt 159999.0
#define ALIGNMENT 32
#define AVX2 1
#define ARCH "AVX2"
#define TOLERANCE 1e-10
#define DELTA 1e-6
#define NDIM3 3
#define DELTA3 1e-4


using Matrix2x2 = std::array<std::array<float, 2>, 2>;
using Matrix3x3 = std::array<std::array<float, 3>, 3>;
using Matrix4x4 = std::array<std::array<float, 4>, 4>;
using MatrixNDIM = std::array<std::array<float, NDIM>, NDIM>;


#include <Tensor.h>
#include <matrix.h>
#include <Metric.h>
#include <Connexion.h>
#include <Grid.h>
#include <GridTensor.h>
#include <Derivatives.h>
#include <Log.h>

typedef struct {
    float x, y, z;
    float lambda;
} GeodesicPoint;

#ifdef  AVX2
    #define VEC_TYPE __m256d
    #define VEC_MUL_PD(a, b) _mm256_mul_pd(a, b)
    #define VEC_SUB_PD(a, b) _mm256_sub_pd(a, b)
    #define VEC_ADD_PD(a, b) _mm256_add_pd(a, b)
    #define VEC_DIV_PD(a, b) _mm256_div_pd(a, b)
    #define VEC_SET1_PD(a)   _mm256_set1_pd(a)
    #define VEC_LOAD_PD(a)   _mm256_load_pd(a)
    #define VEC_STORE_PD(a, b) _mm256_store_pd(a, b)
    #define VEC_CVTSD_F64(a) _mm256_cvtsd_f64(a)
    #define CALCULATE_K_AVX2(k, src_v) \
    for (int mu = 0; mu < 4; mu++) { \
        k[mu] = src_v[mu]; \
        for (int alpha = 0; alpha < 4; alpha++) { \
            __m256d product1 = _mm256_mul_pd(christoffel[mu][alpha][0], \
                                _mm256_mul_pd(src_v[alpha], src_v[0])); \
            __m256d product2 = _mm256_mul_pd(christoffel[mu][alpha][1], \
                                _mm256_mul_pd(src_v[alpha], src_v[1])); \
            __m256d product3 = _mm256_mul_pd(christoffel[mu][alpha][2], \
                                _mm256_mul_pd(src_v[alpha], src_v[2])); \
            __m256d product4 = _mm256_mul_pd(christoffel[mu][alpha][3], \
                                _mm256_mul_pd(src_v[alpha], src_v[3])); \
            k[mu] = _mm256_sub_pd(k[mu], product1); \
            k[mu] = _mm256_sub_pd(k[mu], product2); \
            k[mu] = _mm256_sub_pd(k[mu], product3); \
            k[mu] = _mm256_sub_pd(k[mu], product4); \
        } \
    }

    #define UPDATE_POSITIONS_AVX2(x, v, k, step) \
        for (int mu = 0; mu < 4; mu++) { \
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step, k[mu])); \
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step, k[mu])); \
        }


#else
    #error "AVX2 or AVX512F support required"
#endif


void write_vtk_file(const char *filename);
void store_geodesic_point(float x[4], float lambda);
void geodesic_AVX(__m256d x[4], __m256d v[4], float lambda_max,\
				  __m256d christoffel[4][4][4], __m256d step_size);
void store_geodesic_point_AVX(__m256d x[4], float lambda);

float calculate_impact_parameter(float p_t, float p_phi, float g_tt, float g_tphi, float g_phiphi);
float calculate_emission_angle(float p_r, float p_phi, float g_rr, float g_phiphi);
float b_critique_kerr(float a, int sense);
int  compute_photon_properties(float g[4][4], float p[4]);

/* Problem specific functions */

int Riemann_tensor(const char *metric);
int Geodesics_prob();
int light_geodesics_prob(); 
int Metric_prob();
int grid_setup(); 
void generate_blackhole_image();
void generate_blackhole_shadow();
void evolveADM(Grid::Cell2D& cell, int i, int j, float dt, 
               const std::vector<std::vector<float>>& alpha_grid,
               float r_min, float theta_min, float dr, float dtheta,
               std::ofstream& file, int step);
void calc_gamma_ij_2D(
    int i, int j,
    float r_min, float dr,
    float theta_min, float dtheta,
    Metric &metric_obj,
    Matrix3x3 &gamma3, Matrix3x3 &gamma3_inv);
void calc_gamma_ij(const Vector3& X3D, Matrix3x3& gamma3, Matrix3x3& gamma3_inv);
