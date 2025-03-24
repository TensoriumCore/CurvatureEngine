#include <Geodesics.h>
#include "Spectral.h"
/*
 * Theses functions are used to compute the partial derivative of the metric tensor
 * and the Christoffel symbols, gauge.beta and gauge.alpha gauges conditions
 * All of these use a 2nd, fourth or sixth order finite difference scheme
 *
 * The functions are used in the computation of the ADM/BSSN equations
 *
 * derivatives are computed using the following formulas:
 * \partial_x f = (f_{i+1} - f_{i-1}) / (2 * dx)
 * \partial_y f = (f_{j+1} - f_{j-1}) / (2 * dy)
 * \partial_z f = (f_{k+1} - f_{k-1}) / (2 * dz)
 * \partial_{xx} f = (f_{i+1} - 2 * f_i + f_{i-1}) / dx^2
 * \partial_{yy} f = (f_{j+1} - 2 * f_j + f_{j-1}) / dy^2
 * \partial_{zz} f = (f_{k+1} - 2 * f_k + f_{k-1}) / dz^2
 * \partial_{xy} f = (f_{i+1,j+1} - f_{i+1,j-1} - f_{i-1,j+1} + f_{i-1,j-1}) / (4 * dx * dy)
 * \partial_{xz} f = (f_{i+1,k+1} - f_{i+1,k-1} - f_{i-1,k+1} + f_{i-1,k-1}) / (4 * dx * dz)
 * \partial_{yz} f = (f_{j+1,k+1} - f_{j+1,k-1} - f_{j-1,k+1} + f_{j-1,k-1}) / (4 * dy * dz)
 * \partial_{xxy} f = (f_{i+2} - 2 * f_{i+1} + 2 * f_{i-1} - f_{i-2}) / (12 * dx)
 * \partial_{yyz} f = (f_{j+2} - 2 * f_{j+1} + 2 * f_{j-1} - f_{j-2}) / (12 * dy)
 * \partial_{zzx} f = (f_{k+2} - 2 * f_{k+1} + 2 * f_{k-1} - f_{k-2}) / (12 * dz)
 * \partial_{xxyy} f = (f_{i+1,j+1} - 2 * f_{i+1,j} + f_{i+1,j-1} - f_{i-1,j+1} + 2 * f_{i-1,j} - f_{i-1,j-1}) / (4 * dx * dy)
 * \partial_{yyzz} f = (f_{j+1,k+1} - 2 * f_{j+1,k} + f_{j+1,k-1} - f_{j-1,k+1} + 2 * f_{j-1,k} - f_{j-1,k-1}) / (4 * dy * dz)
 * \partial_{zzxx} f = (f_{k+1,i+1} - 2 * f_{k+1,i} + f_{k+1,i-1} - f_{k-1,i+1} + 2 * f_{k-1,i} - f_{k-1,i-1}) / (4 * dz * dx)
* */


double partialX_alpha(Grid &grid_obj, int i, int j, int k) {
    if (i < 2 || i > NX - 3) {
        return (grid_obj.getCell(i + 1, j, k).gauge.alpha - grid_obj.getCell(i - 1, j, k).gauge.alpha) / (2.0 * DX);
    }

    double alpha_ip2 = grid_obj.getCell(i + 2, j, k).gauge.alpha;
    double alpha_ip1 = grid_obj.getCell(i + 1, j, k).gauge.alpha;
    double alpha_im1 = grid_obj.getCell(i - 1, j, k).gauge.alpha;
    double alpha_im2 = grid_obj.getCell(i - 2, j, k).gauge.alpha;

    return (-alpha_ip2 + 8.0 * alpha_ip1 - 8.0 * alpha_im1 + alpha_im2) / (12.0 * DX);
}


double partialY_alpha(Grid &grid_obj, int i, int j, int k) {
    if (j < 2 || j > NY - 3) {
        return (grid_obj.getCell(i, j + 1, k).gauge.alpha - grid_obj.getCell(i, j - 1, k).gauge.alpha) / (2.0 * DY);
    }
    double alpha_jp2 = grid_obj.getCell(i, j + 2, k).gauge.alpha;
    double alpha_jp1 = grid_obj.getCell(i, j + 1, k).gauge.alpha;
    double alpha_jm1 = grid_obj.getCell(i, j - 1, k).gauge.alpha;
    double alpha_jm2 = grid_obj.getCell(i, j - 2, k).gauge.alpha;

    return (-alpha_jp2 + 8.0 * alpha_jp1 - 8.0 * alpha_jm1 + alpha_jm2) / (12.0 * DY);
}

double partialZ_alpha(Grid &grid_obj, int i, int j, int k) {
    if (k < 2 || k > NZ - 3) {
        return (grid_obj.getCell(i, j, k + 1).gauge.alpha - grid_obj.getCell(i, j, k - 1).gauge.alpha) / (2.0 * DZ);
    }
    double alpha_kp2 = grid_obj.getCell(i, j, k + 2).gauge.alpha;
    double alpha_kp1 = grid_obj.getCell(i, j, k + 1).gauge.alpha;
    double alpha_km1 = grid_obj.getCell(i, j, k - 1).gauge.alpha;
    double alpha_km2 = grid_obj.getCell(i, j, k - 2).gauge.alpha;

    return (-alpha_kp2 + 8.0 * alpha_kp1 - 8.0 * alpha_km1 + alpha_km2) / (12.0 * DZ);
}


double partialXX_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i+1, j, k).gauge.alpha - 2.0 * grid_obj.getCell(i, j, k).gauge.alpha + grid_obj.getCell(i-1, j, k).gauge.alpha) 
           / (DX * DX);
}

double partialYY_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i, j+1, k).gauge.alpha - 2.0 * grid_obj.getCell(i, j, k).gauge.alpha + grid_obj.getCell(i, j-1, k).gauge.alpha) 
           / (DY * DY);
}

double partialZZ_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i, j, k+1).gauge.alpha - 2.0 * grid_obj.getCell(i, j, k).gauge.alpha + grid_obj.getCell(i, j, k-1).gauge.alpha) 
           / (DZ * DZ);
}


double partialX_betacomp(Grid &grid_obj, int i, int j, int k, int comp) {
    return (grid_obj.getCell(i+1, j, k).gauge.beta[comp] - grid_obj.getCell(i-1, j, k).gauge.beta[comp]) / (2.0 * DX);
}

double partialY_betacomp(Grid &grid_obj, int i, int j, int k, int comp) {
    return (grid_obj.getCell(i, j+1, k).gauge.beta[comp] - grid_obj.getCell(i, j-1, k).gauge.beta[comp]) / (2.0 * DY);
}

double partialZ_betacomp(Grid &grid_obj, int i, int j, int k, int comp) {
    return (grid_obj.getCell(i, j, k+1).gauge.beta[comp] - grid_obj.getCell(i, j, k-1).gauge.beta[comp]) / (2.0 * DZ);
}

double partialXY_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i+1, j+1, k).gauge.alpha - grid_obj.getCell(i+1, j-1, k).gauge.alpha
            - grid_obj.getCell(i-1, j+1, k).gauge.alpha + grid_obj.getCell(i-1, j-1, k).gauge.alpha)
           / (4.0 * DX * DY);
}

double partialXZ_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i+1, j, k+1).gauge.alpha - grid_obj.getCell(i+1, j, k-1).gauge.alpha
            - grid_obj.getCell(i-1, j, k+1).gauge.alpha + grid_obj.getCell(i-1, j, k-1).gauge.alpha)
           / (4.0 * DX * DZ);
}

double partialYZ_alpha(Grid &grid_obj, int i, int j, int k) {
    return (grid_obj.getCell(i, j+1, k+1).gauge.alpha - grid_obj.getCell(i, j+1, k-1).gauge.alpha
            - grid_obj.getCell(i, j-1, k+1).gauge.alpha + grid_obj.getCell(i, j-1, k-1).gauge.alpha)
           / (4.0 * DY * DZ);
}

double second_partial_alpha(Grid &grid_obj, int i, int j, int k, int a, int b)
{
	if(a==0 && b==0) return partialXX_alpha(grid_obj, i, j, k);
	if(a==1 && b==1) return partialYY_alpha(grid_obj, i, j, k);
	if(a==2 && b==2) return partialZZ_alpha(grid_obj, i, j, k);

	if((a==0 && b==1)||(a==1 && b==0)) return partialXY_alpha(grid_obj, i, j, k);
	if((a==0 && b==2)||(a==2 && b==0)) return partialXZ_alpha(grid_obj, i, j, k);
	if((a==1 && b==2)||(a==2 && b==1)) return partialYZ_alpha(grid_obj, i, j, k);

	return 0.0;
}


double GridTensor::partialX_gamma(Grid &grid_obj, int i, int j, int k, int a, int b) {
    if (i >= 2 && i <= NX - 3) {
        return fourth_order_diff(
            grid_obj.getCell(i+2, j, k).geom.gamma[a][b],
            grid_obj.getCell(i+1, j, k).geom.gamma[a][b],
            grid_obj.getCell(i-1, j, k).geom.gamma[a][b],
            grid_obj.getCell(i-2, j, k).geom.gamma[a][b],
            DX 
        );
    } else if (i >= 1 && i <= NX - 2) {
        return second_order_diff(
            grid_obj.getCell(i+1, j, k).geom.gamma[a][b],
            grid_obj.getCell(i-1, j, k).geom.gamma[a][b],
            DX 
        );
    } else if (i == 0) {
        return (grid_obj.getCell(i+1, j, k).geom.gamma[a][b] - 
                grid_obj.getCell(i, j, k).geom.gamma[a][b]) / DX;
    } else if (i == NX - 1) {
        return (grid_obj.getCell(i, j, k).geom.gamma[a][b] - 
                grid_obj.getCell(i-1, j, k).geom.gamma[a][b]) / DX;
    }
    return 0.0;
}

double GridTensor::partialY_gamma(Grid &grid_obj, int i, int j, int k, int a, int b) {
    if (j >= 2 && j <= NY - 3) {
        return fourth_order_diff(
            grid_obj.getCell(i, j+2, k).geom.gamma[a][b],
            grid_obj.getCell(i, j+1, k).geom.gamma[a][b],
            grid_obj.getCell(i, j-1, k).geom.gamma[a][b],
            grid_obj.getCell(i, j-2, k).geom.gamma[a][b],
            DY
        );
    } else if (j >= 1 && j <= NY - 2) {
        return second_order_diff(
            grid_obj.getCell(i, j+1, k).geom.gamma[a][b],
            grid_obj.getCell(i, j-1, k).geom.gamma[a][b],
            DY
        );
    } else if (j == 0) {
        return (grid_obj.getCell(i, j+1, k).geom.gamma[a][b] - 
                grid_obj.getCell(i, j, k).geom.gamma[a][b]) / DY;
    } else if (j == NY - 1) {
        return (grid_obj.getCell(i, j, k).geom.gamma[a][b] - 
                grid_obj.getCell(i, j-1, k).geom.gamma[a][b]) / DY;
    }
    return 0.0;
}

/*
 * The partial derivative of the metric tensor with respect to z
 * it uses a 2nd, fourth or sixth order finite difference scheme
 * */

double GridTensor::partialZ_gamma(Grid &grid_obj, int i, int j, int k, int a, int b) {
    if (k >= 2 && k <= NZ - 3) {
        return fourth_order_diff(
            grid_obj.getCell(i, j, k+2).geom.gamma[a][b],
            grid_obj.getCell(i, j, k+1).geom.gamma[a][b],
            grid_obj.getCell(i, j, k-1).geom.gamma[a][b],
            grid_obj.getCell(i, j, k-2).geom.gamma[a][b],
            DZ
        );
    } else if (k >= 1 && k <= NZ - 2) {
        return second_order_diff(
            grid_obj.getCell(i, j, k+1).geom.gamma[a][b],
            grid_obj.getCell(i, j, k-1).geom.gamma[a][b],
            DZ
        );
    } else if (k == 0) {
        return (grid_obj.getCell(i, j, k+1).geom.gamma[a][b] - 
                grid_obj.getCell(i, j, k).geom.gamma[a][b]) / DZ;
    } else if (k == NZ - 1) {
        return (grid_obj.getCell(i, j, k).geom.gamma[a][b] - 
                grid_obj.getCell(i, j, k-1).geom.gamma[a][b]) / DZ;
    }
    return 0.0;
}

/*
 * these functions are used to compute the partial derivative of the extrinsic curvature tensor
 * and the Christoffel symbols, gauge.beta and gauge.alpha gauges conditions
 * All of these use a 2nd, fourth or sixth order finite difference scheme
 * it's working nicely but it's not optimized and probably not the best way to do it
 *
 * I'm currently working on a new method based on the spectral differentiation using Chebyshev polynomials
 * and the barycentric interpolation or the Lagrange interpolation (not sure yet)
 * */


double GridTensor::partialX_Kij(Grid &grid_obj, int i, int j, int k, int a, int b)
{
	/*
	 * The fourth order difference is used when the point is not near the boundary
	 * */
	if (i >= 2 && i <= NX - 3) {
		return fourth_order_diff(
			grid_obj.getCell(i+2, j, k).curv.K[a][b],
			grid_obj.getCell(i+1, j, k).curv.K[a][b],
			grid_obj.getCell(i-1, j, k).curv.K[a][b],
			grid_obj.getCell(i-2, j, k).curv.K[a][b],
			DX 
		);

	/*
	 * The second order difference is used when the point is near the boundary
	 *  This is done to avoid the out of bounds error (boundary conditions are not implemented yet 
	 *  exept for a Cell copy at the boundary, not really useful for now but it's a start)
	 * */

	} else if (i >= 1 && i <= NX - 2) {
		return second_order_diff(
			grid_obj.getCell(i+1, j, k).curv.K[a][b],
			grid_obj.getCell(i-1, j, k).curv.K[a][b],
			DX 
		);
		/*
		 * The first order difference is used when the point is at the boundary
		 * */
	} else if (i == 0) {
		return (grid_obj.getCell(i+1, j, k).curv.K[a][b] - 
				grid_obj.getCell(i, j, k).curv.K[a][b]) / DX;
	} else if (i == NX - 1) {
		return (grid_obj.getCell(i, j, k).curv.K[a][b] - 
				grid_obj.getCell(i-1, j, k).curv.K[a][b]) / DX;
	}
	return 0.0;
}

double GridTensor::partialY_Kij(Grid &grid_obj, int i, int j, int k, int a, int b)
{
	if (j >= 2 && j <= NY - 3) {
		return fourth_order_diff(
			grid_obj.getCell(i, j+2, k).curv.K[a][b],
			grid_obj.getCell(i, j+1, k).curv.K[a][b],
			grid_obj.getCell(i, j-1, k).curv.K[a][b],
			grid_obj.getCell(i, j-2, k).curv.K[a][b],
			DY
		);
	} else if (j >= 1 && j <= NY - 2) {
		return second_order_diff(
			grid_obj.getCell(i, j+1, k).curv.K[a][b],
			grid_obj.getCell(i, j-1, k).curv.K[a][b],
			DY
		);
	} else if (j == 0) {
		return (grid_obj.getCell(i, j+1, k).curv.K[a][b] - 
				grid_obj.getCell(i, j, k).curv.K[a][b]) / DY;
	} else if (j == NY - 1) {
		return (grid_obj.getCell(i, j, k).curv.K[a][b] - 
				grid_obj.getCell(i, j-1, k).curv.K[a][b]) / DY;
	}
	return 0.0;
}


double GridTensor::partialZ_Kij(Grid &grid_obj, int i, int j, int k, int a, int b)
{
	if (k >= 2 && k <= NZ - 3) {
		return fourth_order_diff(
			grid_obj.getCell(i, j, k+2).curv.K[a][b],
			grid_obj.getCell(i, j, k+1).curv.K[a][b],
			grid_obj.getCell(i, j, k-1).curv.K[a][b],
			grid_obj.getCell(i, j, k-2).curv.K[a][b],
			DZ
		);
	} else if (k >= 1 && k <= NZ - 2) {
		return second_order_diff(
			grid_obj.getCell(i, j, k+1).curv.K[a][b],
			grid_obj.getCell(i, j, k-1).curv.K[a][b],
			DZ
		);
	} else if (k == 0) {
		return (grid_obj.getCell(i, j, k+1).curv.K[a][b] - 
				grid_obj.getCell(i, j, k).curv.K[a][b]) / DZ;
	} else if (k == NZ - 1) {
		return (grid_obj.getCell(i, j, k).curv.K[a][b] - 
				grid_obj.getCell(i, j, k-1).curv.K[a][b]) / DZ;
	}
	return 0.0;
}



