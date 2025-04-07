#include <Geodesics.h>

extern float	(*geodesic_points)[5];
extern int		num_points;
int				capacity = 0;
extern float	a;
/* 
	* Write the geodesic points to a VTK file
	* The VTK file will contain the geodesic points and the lambda values
	* The lambda values are the affine parameter values
	* The VTK file can be visualized using Paraview
*/

void write_vtk_file(const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Geodesic Points\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");
    fprintf(file, "POINTS %d float\n", num_points);

    for (int i = 0; i < num_points; ++i)
        fprintf(file, "%f %f %f\n", geodesic_points[i][0], geodesic_points[i][1], geodesic_points[i][2]);

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
    #pragma omp critical
    {
        if (num_points + 16 >= capacity) {
            capacity = (capacity == 0) ? 1000 : capacity * 2;
            float (*new_geodesic_points)[5];
            if (posix_memalign((void **)&new_geodesic_points, ALIGNMENT,
                               capacity * sizeof(*geodesic_points)) != 0)
            {
                fprintf(stderr, "Error: failed to allocate memory for geodesic_points\n");
                exit(EXIT_FAILURE);
            }

            if (geodesic_points) {
                memcpy(new_geodesic_points, geodesic_points,
                       num_points * sizeof(*geodesic_points));
                free(geodesic_points);
            }

            geodesic_points = new_geodesic_points;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < 4; i++) {
        __attribute__((aligned(32))) float r_vals[4], th_vals[4], ph_vals[4];

        VEC_TYPE r     = x[1];
        VEC_TYPE theta = x[2];
        VEC_TYPE phi   = x[3];

        _mm256_store_pd((double *)r_vals,   r);
        _mm256_store_pd((double *)th_vals,  theta);
        _mm256_store_pd((double *)ph_vals,  phi);

        float sin_th[4], cos_th[4], sin_ph[4], cos_ph[4];
        for (int j = 0; j < 4; j++) {
            sin_th[j] = sin(th_vals[j]);
            cos_th[j] = cos(th_vals[j]);
            sin_ph[j] = sin(ph_vals[j]);
            cos_ph[j] = cos(ph_vals[j]);
        }

        for (int j = 0; j < 4; j++) {
            float rr = r_vals[j];
            float xx = rr * sin_th[j] * cos_ph[j];
            float yy = rr * sin_th[j] * sin_ph[j];
            float zz = rr * cos_th[j];

            #pragma omp critical
            {
                geodesic_points[num_points][0] = xx;
                geodesic_points[num_points][1] = yy;
                geodesic_points[num_points][2] = zz;
                geodesic_points[num_points][3] = lambda;
                num_points++;
            }
        }
    }
}
