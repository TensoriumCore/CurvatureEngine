#pragma once

inline constexpr float DX = 0.08f;
inline constexpr float DY = 0.08f;
inline constexpr float DZ = 0.08f;
inline constexpr int NX = 128;
inline constexpr int NY = 128;
inline constexpr int NZ = 128;
inline constexpr int GHOST = 2;
inline constexpr int NX_TOTAL = NX + 2 * GHOST;
inline constexpr int NY_TOTAL = NY + 2 * GHOST;
inline constexpr int NZ_TOTAL = NZ + 2 * GHOST;
