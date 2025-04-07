#include <Geodesics.h>

float Matrix::determinant2x2(const std::array<std::array<float, 2>, 2>& mat) {
    return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

float Matrix::determinant3x3(const Matrix3x3& mat) {
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
         - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
         + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

float Matrix::determinant4x4(const Matrix4x4& mat) {
	Matrix3x3 minor{};
    float det = 0.0;
    for (int i = 0; i < NDIM; i++) {
        int subi = 0; 
        for (int j = 1; j < NDIM; j++) {
            int subj = 0;
            for (int k = 0; k < NDIM; k++) {
                if (k == i)
                    continue;
                minor[subi][subj] = mat[j][k];
                subj++;
            }
            subi++;
        }
        det += pow(-1, i) * mat[0][i] * determinant3x3(minor);
    }
    return det;
}

void Matrix::cofactor(const Matrix4x4& mat, Matrix4x4& cofactorMat) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Matrix3x3 minor{};
            int subi = 0;
            for (int x = 0; x < 4; x++) {
                if (x == i) continue;
                int subj = 0;
                for (int y = 0; y < 4; y++) {
                    if (y == j) continue;
                    minor[subi][subj] = mat[x][y];
                    subj++;
                }
                subi++;
            }
            cofactorMat[i][j] = (i + j) % 2 == 0 ? determinant3x3(minor) : -determinant3x3(minor);
        }
    }
}

void Matrix::cofactor3x3(const Matrix3x3& mat, Matrix3x3& cofactorMat) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Matrix2x2 minor{};
            int subi = 0;
            for (int x = 0; x < 3; x++) {
                if (x == i) continue; 
                int subj = 0;
                for (int y = 0; y < 3; y++) {
                    if (y == j) continue; 
                    minor[subi][subj] = mat[x][y];
                    subj++;
                }
                subi++;
            }
            cofactorMat[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * determinant2x2(minor);
        }
    }
}

void Matrix::transpose(const MatrixNDIM& mat, MatrixNDIM& transposed) {
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            transposed[j][i] = mat[i][j];
        }
    }
}

void Matrix::transpose3x3(const Matrix3x3& mat, Matrix3x3& transposed) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            transposed[j][i] = mat[i][j];
        }
    }
}

int Matrix::inverse_matrix(const MatrixNDIM& mat, MatrixNDIM& inverse) {
    float det = determinant4x4(mat);
    if (fabs(det) < 1e-10) {
        return 0; 
    }

    MatrixNDIM cofactorMat{}, adjugate{};
    cofactor(mat, cofactorMat);
    transpose(cofactorMat, adjugate);

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            inverse[i][j] = adjugate[i][j] / det;
        }
    }
    return 1;
}

int Matrix::inverse_3x3(const Matrix3x3& mat, Matrix3x3& inv) {
    float det = determinant3x3(mat);
    if (fabs(det) < 1e-6) {
        return 0; 
    }

    Matrix3x3 cofactorMat{}, adjugate{};
    cofactor3x3(mat, cofactorMat);
    transpose3x3(cofactorMat, adjugate);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inv[i][j] = adjugate[i][j] / det;
        }
    }
    return 1;
}

void Matrix::check_inverse_3x3(const Matrix3x3& mat, const Matrix3x3& inv) {
    Matrix3x3 product{};
    constexpr Matrix3x3 identity = {};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            product[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                product[i][j] += mat[i][k] * inv[k][j];
            }
        }
    }

    constexpr float TOL = 1e-10;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (fabs(product[i][j] - identity[i][j]) > TOL) {
                std::cerr << "Inversion 3x3 check failed at (" << i << ", " << j << "): "
                          << product[i][j] << " != " << identity[i][j] << std::endl;
            }
        }
    }
}



void Matrix::check_inverse(const MatrixNDIM& gcov, const MatrixNDIM& g_inv) {
    MatrixNDIM identity{};
    for (int i = 0; i < NDIM; i++) {
        identity[i][i] = 1.0;
    }

    MatrixNDIM product{};
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            for (int k = 0; k < NDIM; k++) {
                product[i][j] += gcov[i][k] * g_inv[k][j];
            }
        }
    }

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            if (std::fabs(product[i][j] - identity[i][j]) > TOLERANCE) {
                std::cerr << "Matrix inversion check failed at element (" << i << ", " << j << "): " 
                          << product[i][j] << std::endl;
            }
        }
    }
	print_matrix("identity", identity);
}


void Matrix::print_matrix(const char* name, const MatrixNDIM& mat) {
    std::cout << name << " =\n";
    for (const auto& row : mat) {
        for (const auto& elem : row) {
            std::cout << std::fixed << std::setprecision(12) << elem << " ";
        }
        std::cout << std::endl;
    }
}

void Matrix::print_matrix_3x3(const char* name, const Matrix3x3& mat) {
    std::cout << name << " =\n";
    for (const auto& row : mat) {
        for (const auto& elem : row) {
            std::cout << std::fixed << std::setprecision(12) << elem << " ";
        }
        std::cout << std::endl;
    }
}
