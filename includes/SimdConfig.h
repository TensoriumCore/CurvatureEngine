#pragma once

#include <cstddef>
#include <cstdint>

#ifndef CURVATUREENGINE_ENABLE_AVX2
#define CURVATUREENGINE_ENABLE_AVX2 0
#endif

#ifndef CURVATUREENGINE_ENABLE_NEON
#define CURVATUREENGINE_ENABLE_NEON 0
#endif

#if CURVATUREENGINE_ENABLE_AVX2 && defined(__x86_64__) && defined(__AVX2__)
#define CURVATUREENGINE_TARGET_AVX2 1
#else
#define CURVATUREENGINE_TARGET_AVX2 0
#endif

#if CURVATUREENGINE_ENABLE_NEON && defined(__aarch64__)
#define CURVATUREENGINE_TARGET_NEON 1
#else
#define CURVATUREENGINE_TARGET_NEON 0
#endif

namespace curvatureengine::simd {

#if CURVATUREENGINE_TARGET_AVX2
#include <immintrin.h>
using Vec4d = __m256d;

inline Vec4d broadcast(double value) { return _mm256_set1_pd(value); }
inline Vec4d zero() { return _mm256_setzero_pd(); }
inline Vec4d load(const double *ptr) { return _mm256_load_pd(ptr); }
inline void store(double *ptr, Vec4d value) { _mm256_store_pd(ptr, value); }
inline Vec4d add(Vec4d a, Vec4d b) { return _mm256_add_pd(a, b); }
inline Vec4d sub(Vec4d a, Vec4d b) { return _mm256_sub_pd(a, b); }
inline Vec4d mul(Vec4d a, Vec4d b) { return _mm256_mul_pd(a, b); }
inline Vec4d div(Vec4d a, Vec4d b) { return _mm256_div_pd(a, b); }
inline Vec4d fmadd(Vec4d a, Vec4d b, Vec4d c) {
  return _mm256_fmadd_pd(a, b, c);
}
inline Vec4d fnmadd(Vec4d a, Vec4d b, Vec4d c) {
  return _mm256_fnmadd_pd(a, b, c);
}
inline double lane0(Vec4d value) { return _mm256_cvtsd_f64(value); }
#elif CURVATUREENGINE_TARGET_NEON
#include <arm_neon.h>
struct alignas(32) Vec4d {
  float64x2_t lo;
  float64x2_t hi;
};

inline Vec4d broadcast(double value) {
  return {vdupq_n_f64(value), vdupq_n_f64(value)};
}

inline Vec4d zero() { return broadcast(0.0); }

inline Vec4d load(const double *ptr) {
  return {vld1q_f64(ptr), vld1q_f64(ptr + 2)};
}

inline void store(double *ptr, const Vec4d &value) {
  vst1q_f64(ptr, value.lo);
  vst1q_f64(ptr + 2, value.hi);
}

inline Vec4d add(const Vec4d &a, const Vec4d &b) {
  return {vaddq_f64(a.lo, b.lo), vaddq_f64(a.hi, b.hi)};
}

inline Vec4d sub(const Vec4d &a, const Vec4d &b) {
  return {vsubq_f64(a.lo, b.lo), vsubq_f64(a.hi, b.hi)};
}

inline Vec4d mul(const Vec4d &a, const Vec4d &b) {
  return {vmulq_f64(a.lo, b.lo), vmulq_f64(a.hi, b.hi)};
}

inline Vec4d div(const Vec4d &a, const Vec4d &b) {
  return {vdivq_f64(a.lo, b.lo), vdivq_f64(a.hi, b.hi)};
}

inline Vec4d fmadd(const Vec4d &a, const Vec4d &b, const Vec4d &c) {
#if defined(__ARM_FEATURE_FMA)
  return {vfmaq_f64(c.lo, a.lo, b.lo), vfmaq_f64(c.hi, a.hi, b.hi)};
#else
  return add(mul(a, b), c);
#endif
}

inline Vec4d fnmadd(const Vec4d &a, const Vec4d &b, const Vec4d &c) {
  return sub(c, mul(a, b));
}

inline double lane0(const Vec4d &value) { return vgetq_lane_f64(value.lo, 0); }
#else
struct alignas(32) Vec4d {
  double lane[4];
};

inline Vec4d broadcast(double value) {
  return Vec4d{{value, value, value, value}};
}
inline Vec4d zero() { return broadcast(0.0); }
inline Vec4d load(const double *ptr) {
  return Vec4d{{ptr[0], ptr[1], ptr[2], ptr[3]}};
}
inline void store(double *ptr, const Vec4d &value) {
  for (int i = 0; i < 4; ++i) {
    ptr[i] = value.lane[i];
  }
}
inline Vec4d add(const Vec4d &a, const Vec4d &b) {
  return Vec4d{{a.lane[0] + b.lane[0], a.lane[1] + b.lane[1],
                a.lane[2] + b.lane[2], a.lane[3] + b.lane[3]}};
}
inline Vec4d sub(const Vec4d &a, const Vec4d &b) {
  return Vec4d{{a.lane[0] - b.lane[0], a.lane[1] - b.lane[1],
                a.lane[2] - b.lane[2], a.lane[3] - b.lane[3]}};
}
inline Vec4d mul(const Vec4d &a, const Vec4d &b) {
  return Vec4d{{a.lane[0] * b.lane[0], a.lane[1] * b.lane[1],
                a.lane[2] * b.lane[2], a.lane[3] * b.lane[3]}};
}
inline Vec4d div(const Vec4d &a, const Vec4d &b) {
  return Vec4d{{a.lane[0] / b.lane[0], a.lane[1] / b.lane[1],
                a.lane[2] / b.lane[2], a.lane[3] / b.lane[3]}};
}
inline Vec4d fmadd(const Vec4d &a, const Vec4d &b, const Vec4d &c) {
  return add(mul(a, b), c);
}
inline Vec4d fnmadd(const Vec4d &a, const Vec4d &b, const Vec4d &c) {
  return sub(c, mul(a, b));
}
inline double lane0(const Vec4d &value) { return value.lane[0]; }
#endif

} // namespace curvatureengine::simd
