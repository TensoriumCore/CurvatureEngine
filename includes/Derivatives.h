#pragma once
#include <Grid.h> 
class Grid;

template <class T>
inline double fourth_order_diff(const T &plus2, const T &plus1, 
                         const T &minus1, const T &minus2, double dx)
{
    return (-plus2 + 8.0*plus1 - 8.0*minus1 + minus2) / (12.0*dx);
}


template <class T>
inline double second_order_diff(const T &plus1, const T &minus1, double dx)
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

        double partialX_alpha(int i, int j, int k);
        double partialY_alpha(int i, int j, int k);
        double partialZ_alpha(int i, int j, int k);

        double partialXX_alpha(int i, int j, int k);
        double partialYY_alpha(int i, int j, int k);
        double partialZZ_alpha(int i, int j, int k);

        double partialX_betacomp(int i, int j, int k, int comp);
        double partialY_betacomp(int i, int j, int k, int comp);
        double partialZ_betacomp(int i, int j, int k, int comp);

        double partialXY_alpha(int i, int j, int k);
        double partialXZ_alpha(int i, int j, int k);
        double partialYZ_alpha(int i, int j, int k);

        double second_partial_alpha(int i, int j, int k, int a, int b);
};

template <typename Func>
double partialX(Grid &grid_obj, int i, int j, int k, Func f) {
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
double partialY(Grid &grid_obj, int i, int j, int k, Func f) {
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
double partialZ(Grid &grid_obj, int i, int j, int k, Func f) {
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
double partial_m(Grid &grid_obj, int i, int j, int k, int dim, Getter getter) {
    if (dim == 0) return partialX(grid_obj, i, j, k, getter);
    if (dim == 1) return partialY(grid_obj, i, j, k, getter);
    if (dim == 2) return partialZ(grid_obj, i, j, k, getter);
    return 0.0;
}


