#pragma once

#include "Geodesics.h"

class Matrix {
	public:
		float determinant3x3(const Matrix3x3& mat);
		float determinant4x4(const Matrix4x4& mat);
		void cofactor(const MatrixNDIM& mat, MatrixNDIM& cofactorMat);
		void transpose(const MatrixNDIM& mat, MatrixNDIM& transposed);
		int inverse_matrix(const MatrixNDIM& mat, MatrixNDIM& inverse);
		void check_inverse(const MatrixNDIM& gcov, const MatrixNDIM& gcon);
		void print_matrix(const char *name, const MatrixNDIM& matrix);
		int inverse_3x3(const Matrix3x3& mat, Matrix3x3& inv);
		void check_inverse_3x3(const Matrix3x3& mat, const Matrix3x3& inv);
		void print_matrix_3x3(const char* name, const Matrix3x3& mat);
		void transpose3x3(const Matrix3x3& mat, Matrix3x3& transposed);
		void cofactor3x3(const Matrix3x3& mat, Matrix3x3& cofactorMat);
		void check_cofactor3x3(const Matrix3x3& mat, const Matrix3x3& cofactorMat);
		float determinant2x2(const std::array<std::array<float, 2>, 2>& mat);
};
