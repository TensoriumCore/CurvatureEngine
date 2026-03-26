
#pragma once

#include "core/Constants.h"
#include "core/GeodesicIntegrator.h"
#include "core/Types.h"

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

#include "Connexion.h"
#include "Derivatives.h"
#include "Grid.h"
#include "GridTensor.h"
#include "Log.h"
#include "Metric.h"
#include "Tensor.h"
#include "matrix.h"

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
int shadow_prob(); 
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
