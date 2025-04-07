#include <Geodesics.h>

extern float a;

void geodesic_AVX(__m256d x[4], __m256d v[4], float lambda_max,
                   __m256d christoffel[4][4][4], __m256d step_size) 
{
    __attribute__((aligned(32))) __m256d k1_x[4], k1_v[4];
    __attribute__((aligned(32))) __m256d k2_x[4], k2_v[4];
    __attribute__((aligned(32))) __m256d k3_x[4], k3_v[4];
    __attribute__((aligned(32))) __m256d k4_x[4], k4_v[4];
    __attribute__((aligned(32))) __m256d temp_x[4], temp_v[4];
    __attribute__((aligned(32)))  float lambda = 0.0;
	float min_r = 1.0 + sqrt(1.0 - a * a) + 1;
    while (lambda < lambda_max) {

        for (int mu = 0; mu < 4; mu++) {
            k1_x[mu] = v[mu];
            k1_v[mu] = _mm256_setzero_pd();
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d prod = _mm256_mul_pd(v[alpha], v[beta]);
                    k1_v[mu] = _mm256_fnmadd_pd(christoffel[mu][alpha][beta], prod, k1_v[mu]);
                }
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            __m256d half_step = _mm256_mul_pd(step_size, _mm256_set1_pd(0.5));
            temp_x[mu] = _mm256_fmadd_pd(half_step, k1_x[mu], x[mu]);
            temp_v[mu] = _mm256_fmadd_pd(half_step, k1_v[mu], v[mu]);
        }

        for (int mu = 0; mu < 4; mu++) {
            k2_x[mu] = temp_v[mu];
            k2_v[mu] = _mm256_setzero_pd();
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d prod = _mm256_mul_pd(temp_v[alpha], temp_v[beta]);
                    k2_v[mu] = _mm256_fnmadd_pd(christoffel[mu][alpha][beta], prod, k2_v[mu]);
                }
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            __m256d half_step = _mm256_mul_pd(step_size, _mm256_set1_pd(0.5));
            temp_x[mu] = _mm256_fmadd_pd(half_step, k2_x[mu], x[mu]);
            temp_v[mu] = _mm256_fmadd_pd(half_step, k2_v[mu], v[mu]);
        }

        for (int mu = 0; mu < 4; mu++) {
            k3_x[mu] = temp_v[mu];
            k3_v[mu] = _mm256_setzero_pd();
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d prod = _mm256_mul_pd(temp_v[alpha], temp_v[beta]);
                    k3_v[mu] = _mm256_fnmadd_pd(christoffel[mu][alpha][beta], prod, k3_v[mu]);
                }
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            temp_x[mu] = _mm256_fmadd_pd(step_size, k3_x[mu], x[mu]);
            temp_v[mu] = _mm256_fmadd_pd(step_size, k3_v[mu], v[mu]);
        }

        for (int mu = 0; mu < 4; mu++) {
            k4_x[mu] = temp_v[mu];
            k4_v[mu] = _mm256_setzero_pd();
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d prod = _mm256_mul_pd(temp_v[alpha], temp_v[beta]);
                    k4_v[mu] = _mm256_fnmadd_pd(christoffel[mu][alpha][beta], prod, k4_v[mu]);
                }
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            __m256d two = _mm256_set1_pd(2.0);
            __m256d six = _mm256_set1_pd(6.0);
            __m256d sum_x = _mm256_add_pd(k1_x[mu],
                               _mm256_add_pd(_mm256_mul_pd(two, k2_x[mu]),
                               _mm256_add_pd(_mm256_mul_pd(two, k3_x[mu]), k4_x[mu])));
            __m256d sum_v = _mm256_add_pd(k1_v[mu],
                               _mm256_add_pd(_mm256_mul_pd(two, k2_v[mu]),
                               _mm256_add_pd(_mm256_mul_pd(two, k3_v[mu]), k4_v[mu])));
            
            x[mu] = _mm256_fmadd_pd(step_size, _mm256_div_pd(sum_x, six), x[mu]);
            v[mu] = _mm256_fmadd_pd(step_size, _mm256_div_pd(sum_v, six), v[mu]);
        }
		if (_mm256_cvtsd_f64(x[1]) < min_r) {	
			printf("particle has crossed the horizon at r = %f\n", _mm256_cvtsd_f64(x[1]));
		}
        lambda += _mm256_cvtsd_f64(step_size);
        store_geodesic_point_AVX(x, lambda);
    }
}
