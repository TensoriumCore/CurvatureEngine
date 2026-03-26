#pragma once

#include "SimdConfig.h"

inline constexpr float C = 1.0f;
inline constexpr float G = 1.0f;
inline constexpr float M = 1.0f;
inline constexpr int BLOCK_SIZE = 1024;
inline constexpr int BUFFER_SIZE = 1024;
inline constexpr float SMALL = 1.0e-40f;
inline constexpr int NDIM = 4;
inline constexpr int DIM3 = 3;
inline constexpr int NDIM3 = DIM3;
inline constexpr float DT = 0.0000005f;
inline constexpr float max_dt = 5999.0f;
inline constexpr int ALIGNMENT = 32;
inline constexpr double TOLERANCE = 1.0e-10;
inline constexpr float DELTA = 1.0e-6f;
inline constexpr float DELTA3 = 1.0e-4f;

#if CURVATUREENGINE_TARGET_AVX2
inline constexpr char ARCH[] = "AVX2";
#elif CURVATUREENGINE_TARGET_NEON
inline constexpr char ARCH[] = "NEON";
#else
inline constexpr char ARCH[] = "SCALAR";
#endif
