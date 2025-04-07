#pragma once
#include <Grid.h> 
class Grid;


inline float safe_dx(float dx) {
    return (dx > 1e-12) ? dx : 1e-12;  // Évite les divisions instables
}



template <class T>
__attribute__((always_inline))
inline float fourth_order_diff(const T &plus2, const T &plus1, 
                                const T &minus1, const T &minus2, float dx)
{
    dx = safe_dx(dx);
    float t1 = std::fma(8.0, plus1, -plus2);   
    float t2 = std::fma(-8.0, minus1, minus2); 
    float numerator = t1 + t2;
    return numerator / (12.0 * dx);
}


template <class T>
__attribute__((always_inline))
inline float second_order_diff(const T &plus1, const T &minus1, float dx)
{
    return (plus1 - minus1)/(2.0*dx);
}


class Derivatives
{
    public:
        Derivatives() = default;
        ~Derivatives() = default;
        friend class Grid;
        friend class GridTensor;
        Grid* grid;

        float partialX_alpha(int i, int j, int k);
        float partialY_alpha(int i, int j, int k);
        float partialZ_alpha(int i, int j, int k);

        float partialXX_alpha(int i, int j, int k);
        float partialYY_alpha(int i, int j, int k);
        float partialZZ_alpha(int i, int j, int k);

        float partialX_betacomp(int i, int j, int k, int comp);
        float partialY_betacomp(int i, int j, int k, int comp);
        float partialZ_betacomp(int i, int j, int k, int comp);

        float partialXY_alpha(int i, int j, int k);
        float partialXZ_alpha(int i, int j, int k);
        float partialYZ_alpha(int i, int j, int k);

		float second_partial_alpha(int i, int j, int k, int a, int b);
};


void spectral_derivative_3D(const std::vector<float> &f_in,
				std::vector<float> &f_out,
				int N_x, int N_y, int N_z,
				float dx, float dy, float dz,
				int dim);


template <typename Func>
float partialX(Grid &grid_obj, int i, int j, int k, Func f) {
	if (i < 0 || i >= NX) return 0.0;
	if (j < 0 || j >= NY) return 0.0;
	if (k < 0 || k >= NZ) return 0.0;
	if (i >= 2 && i <= NX - 3) {
		return fourth_order_diff(
				f(grid_obj.getCell(i+2, j, k)),
				f(grid_obj.getCell(i+1, j, k)),
				f(grid_obj.getCell(i-1, j, k)),
				f(grid_obj.getCell(i-2, j, k)),
				DX
				);
	} else if (i >= 1 && i <= NX - 2) {
		return second_order_diff(
				f(grid_obj.getCell(i+1, j, k)),
				f(grid_obj.getCell(i-1, j, k)),
				DX
				);
	} else if (i == 0) {
		return (f(grid_obj.getCell(i+1, j, k)) - f(grid_obj.getCell(i, j, k))) / DX;
	} else if (i == NX - 1) {
		return (f(grid_obj.getCell(i, j, k)) - f(grid_obj.getCell(i-1, j, k))) / DX;
	}
	return 0.0;
}


template <typename Func>
float partialY(Grid &grid_obj, int i, int j, int k, Func f) {
	if (j < 0 || j >= NY) return 0.0;
	if (i < 0 || i >= NX) return 0.0;
	if (k < 0 || k >= NZ) return 0.0;
	if (j >= 2 && j <= NY - 3) {
		return fourth_order_diff(
				f(grid_obj.getCell(i, j+2, k)),
				f(grid_obj.getCell(i, j+1, k)),
			f(grid_obj.getCell(i, j-1, k)),
			f(grid_obj.getCell(i, j-2, k)),
			DY
		);
	} else if (j >= 1 && j <= NY - 2) {
		return second_order_diff(
			f(grid_obj.getCell(i, j+1, k)),
			f(grid_obj.getCell(i, j-1, k)),
			DY
		);
	} else if (j == 0) {
		return (f(grid_obj.getCell(i, j+1, k)) - f(grid_obj.getCell(i, j, k))) / DY;
	} else if (j == NY - 1) {
		return (f(grid_obj.getCell(i, j, k)) - f(grid_obj.getCell(i, j-1, k))) / DY;
	}
	return 0.0;
}

template <typename Func>
float partialZ(Grid &grid_obj, int i, int j, int k, Func f) {
	if (k < 0 || k >= NZ) return 0.0;
	if (i < 0 || i >= NX) return 0.0;
	if (j < 0 || j >= NY) return 0.0;
	if (k >= 2 && k <= NZ - 3) {
		return fourth_order_diff(
			f(grid_obj.getCell(i, j, k+2)),
			f(grid_obj.getCell(i, j, k+1)),
			f(grid_obj.getCell(i, j, k-1)),
			f(grid_obj.getCell(i, j, k-2)),
			DZ
		);
	} else if (k >= 1 && k <= NZ - 2) {
		return second_order_diff(
			f(grid_obj.getCell(i, j, k+1)),
			f(grid_obj.getCell(i, j, k-1)),
			DZ
		);
	} else if (k == 0) {
		return (f(grid_obj.getCell(i, j, k+1)) - f(grid_obj.getCell(i, j, k))) / DZ;
	} else if (k == NZ - 1) {
		return (f(grid_obj.getCell(i, j, k)) - f(grid_obj.getCell(i, j, k-1))) / DZ;
	}
	return 0.0;
}


template <typename Getter>
float partial_m(Grid &grid_obj, int i, int j, int k, int dim, Getter getter) {
    if (dim == 0) return partialX(grid_obj, i, j, k, getter);
    if (dim == 1) return partialY(grid_obj, i, j, k, getter);
    if (dim == 2) return partialZ(grid_obj, i, j, k, getter);
    return 0.0;
}

template <typename Getter>
float second_partial(Grid &grid_obj, int i, int j, int k, int a, int b, Getter getter) {
    if (a == b) {
        if (a == 0 && i >= 1 && i <= NX - 2) {
            return (getter(grid_obj.getCell(i+1, j, k)) - 2.0 * getter(grid_obj.getCell(i, j, k)) + getter(grid_obj.getCell(i-1, j, k))) / (DX * DX);
        } else if (a == 1 && j >= 1 && j <= NY - 2) {
            return (getter(grid_obj.getCell(i, j+1, k)) - 2.0 * getter(grid_obj.getCell(i, j, k)) + getter(grid_obj.getCell(i, j-1, k))) / (DY * DY);
        } else if (a == 2 && k >= 1 && k <= NZ - 2) {
            return (getter(grid_obj.getCell(i, j, k+1)) - 2.0 * getter(grid_obj.getCell(i, j, k)) + getter(grid_obj.getCell(i, j, k-1))) / (DZ * DZ);
        }
    } else {
        int ip = i, im = i, jp = j, jm = j, kp = k, km = k;
        if (a == 0) { ip = i + 1; im = i - 1; }
        if (a == 1) { jp = j + 1; jm = j - 1; }
        if (a == 2) { kp = k + 1; km = k - 1; }

        if (b == 0) { ip = i + 1; im = i - 1; }
        if (b == 1) { jp = j + 1; jm = j - 1; }
        if (b == 2) { kp = k + 1; km = k - 1; }

        float dx_a = (a == 0) ? DX : (a == 1) ? DY : DZ;
        float dx_b = (b == 0) ? DX : (b == 1) ? DY : DZ;

        dx_a = safe_dx(dx_a);
        dx_b = safe_dx(dx_b);

        return (getter(grid_obj.getCell(ip, jp, kp)) - getter(grid_obj.getCell(ip, jm, km))
              - getter(grid_obj.getCell(im, jp, kp)) + getter(grid_obj.getCell(im, jm, km)))
              / (4.0 * dx_a * dx_b);
    }

    return 0.0;
}

