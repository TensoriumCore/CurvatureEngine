#pragma once

#include "bssn/GridConfig.h"
#include "core/Types.h"

#include <vector>

struct alignas(32) Geometry {
    Matrix3x3 gamma;
    Matrix3x3 gamma_inv;
    Matrix3x3 Ricci;
    Matrix3x3 gamma0;
    float tilde_gamma[3][3];
    float tilde_gamma0[3][3];
    float tildgamma_inv[3][3];
    float dt_tilde_gamma[3][3];
    float dt_tildeGamma[3];
};

struct alignas(32) Connection {
    Tensor3D Gamma3;
    float tildeGamma[3];
    float Christoffel[3][3][3];
};

struct alignas(32) ExtrinsicCurvature {
    Matrix3x3 K;
    Matrix3x3 K0;
    float K_trace;
    float dt_K_trace;
    float H;
    float dKt[3][3];
};

struct alignas(32) AtildeVars {
    float Atilde[3][3];
    float Atilde0[3][3];
    float dt_Atilde[3][3];
    float AtildeStage[4][3][3];
};

struct alignas(32) Gauge {
    float alpha;
    float alpha0;
    float alphaStage[4];
    Vector3 dalpha_dx;

    float beta[3];
    float beta0[3];
    float betaStage[4][3];
    float Bstage[4][3];
    float dt_beta[3];
    float dt_alpha;
};

struct Matter {
    float rho;
    float momentum[3];
    float hamiltonian;
    float p;
    float vx, vy, vz;
    float T[4][4];
};

struct alignas(32) BSSNCell {
    Geometry geom;
    Connection conn;
    ExtrinsicCurvature curv;
    AtildeVars atilde;
    Gauge gauge;
    Matter matter;
    float t;
    float ADMmass;
    float chi;
    float dt_chi;
    float gammaStage[4][3][3];
    float KStage[4][3][3];
    float dgt[3][3];
};

struct alignas(32) BSSNGridStorage {
    BSSNCell* cells;

    float* alpha;
    float* chi;
    float* Atilde[3][3];

    int nx, ny, nz;

    inline int idx(int i, int j, int k) const {
        return i * ny * nz + j * nz + k;
    }
};

using BSSNCellGrid =
    std::vector<std::vector<std::vector<BSSNCell>>>;
