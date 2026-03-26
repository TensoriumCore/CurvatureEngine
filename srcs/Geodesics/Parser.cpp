#include "core/GeodesicIntegrator.h"

#include <cstdio>
#include <atomic>
#include <cmath>

#include "app/RuntimeState.h"

std::atomic<int> global_idx{0}; 
const int MAX_SAFE_CAPACITY = 10000000; 
int capacity = 0;
/*
 * Write the geodesic points to a VTK file
 * The VTK file will contain the geodesic points and the lambda values
 * The lambda values are the affine parameter values
 * The VTK file can be visualized using Paraview
 */

void write_vtk_file(const char *filename) {
  FILE *file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", filename);
    return;
  }

  fprintf(file, "# vtk DataFile Version 3.0\n");
  fprintf(file, "Geodesic Points\n");
  fprintf(file, "ASCII\n");
  fprintf(file, "DATASET POLYDATA\n");
  fprintf(file, "POINTS %d float\n", num_points);

  for (int i = 0; i < num_points; ++i)
    fprintf(file, "%f %f %f\n", geodesic_points[i][0], geodesic_points[i][1],
            geodesic_points[i][2]);

  if (num_points > 1) {
    fprintf(file, "LINES %d %d\n", num_points - 1, 3 * (num_points - 1));
    for (int i = 0; i < num_points - 1; ++i)
      fprintf(file, "2 %d %d\n", i, i + 1);
  }

  fprintf(file, "POINT_DATA %d\n", num_points);
  fprintf(file, "SCALARS lambda float\n");
  fprintf(file, "LOOKUP_TABLE default\n");

  for (int i = 0; i < num_points; ++i)
    fprintf(file, "%f\n", geodesic_points[i][3]);

  printf("Number of points: %d\n", num_points);
  printf("VTK file %s has been written\n", filename);

  fclose(file);
}



void store_geodesic_point_AVX(VEC_TYPE x[4], float lambda) {
  alignas(32) double r_vals[4], th_vals[4], ph_vals[4];
  curvatureengine::simd::store(r_vals, x[1]);
  curvatureengine::simd::store(th_vals, x[2]);
  curvatureengine::simd::store(ph_vals, x[3]);

  double sin_th[4], cos_th[4], sin_ph[4], cos_ph[4];
  for (int j = 0; j < 4; j++) {
    sin_th[j] = std::sin(th_vals[j]);
    cos_th[j] = std::cos(th_vals[j]);
    sin_ph[j] = std::sin(ph_vals[j]);
    cos_ph[j] = std::cos(ph_vals[j]);
  }

  int idx = global_idx.fetch_add(4, std::memory_order_relaxed);

  if (idx + 4 >= MAX_SAFE_CAPACITY) return;

  for (int j = 0; j < 4; j++) {
    double rr = r_vals[j];
    double xx = rr * sin_th[j] * cos_ph[j];
    double yy = rr * sin_th[j] * sin_ph[j];
    double zz = rr * cos_th[j];

    int current = idx + j;
    geodesic_points[current][0] = static_cast<float>(xx);
    geodesic_points[current][1] = static_cast<float>(yy);
    geodesic_points[current][2] = static_cast<float>(zz);
    geodesic_points[current][3] = lambda;
  }
}
