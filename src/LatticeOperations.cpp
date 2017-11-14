//
// Created by Ilya Bryukhanov on 01.11.2017.
//

#include "LatticeOperations.h"

namespace NTBuilder {

    crystal_dir operator*(const crystal_dir& dir, double num) {
        crystal_dir res = dir;
        for (int i = 0; i < 3; ++i) {
            res[i] *= num;
        }
        return res;
    }

    crystal_dir operator+(const crystal_dir& dir1, const crystal_dir& dir2) {
        crystal_dir res = dir1;
        for (int i = 0; i < 3; ++i) {
            res[i] += dir2[i];
        }
        return res;
    }

    crystal_dir operator-(const crystal_dir& dir1, const crystal_dir& dir2) {
        crystal_dir res = dir1;
        for (int i = 0; i < 3; ++i) {
            res[i] -= dir2[i];
        }
        return res;
    }

    int Scalar(r_crystal_dir x, r_crystal_dir y) {
        return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
    }

    bool CheckPairwiseOrthogonality(r_crystal_dir x, r_crystal_dir y, r_crystal_dir z) {
        return !(Scalar(x, y) || Scalar(y, z) || Scalar(x, z));
    }

    bool CheckRightHanded(r_crystal_dir x, r_crystal_dir y, r_crystal_dir z) {
        int xy0 = x[1] * y[2] - x[2] * y[1];
        int xy1 = x[2] * y[0] - x[0] * y[2];
        int xy2 = x[0] * y[1] - x[1] * y[0];
        return xy0 * z[0] + xy1 * z[1] + xy2 * z[2] > 0;
    }

    double getLength(r_crystal_dir x) {
        return sqrt(static_cast<double>(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
    }

    matrix CreateZeroMatrix3x3() {
        matrix a(3);
        for (int i = 0; i < 3; ++i) {
            a[i].resize(3);
        }
        return a;
    }

    matrix TransposeMatrix3x3(const matrix& m) {
        matrix a = CreateZeroMatrix3x3();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                a[i][j] = m[j][i];
            }
        }
        return a;
    }
    matrix Norm3VectorsRows(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal) {
        matrix m(3);
        for (int i = 0; i < 3; ++i) {
            m[i].resize(3);
        }

        double x_length = getLength(x);
        double y_length = getLength(y);
        double normal_length = getLength(normal);
        for (int i = 0; i < 3; ++i) {
            m[0][i] = x[i] / x_length;
            m[1][i] = y[i] / y_length;
            m[2][i] = normal[i] / normal_length;
        }
        return m;
    }

    matrix Norm3VectorsCols(r_crystal_dir x, r_crystal_dir y, r_crystal_dir normal) {
        matrix m(3);
        for (int i = 0; i < 3; ++i) {
            m[i].resize(3);
        }

        double x_length = getLength(x);
        double y_length = getLength(y);
        double normal_length = getLength(normal);
        for (int i = 0; i < 3; ++i) {
            m[i][0] = x[i] / x_length;
            m[i][1] = y[i] / y_length;
            m[i][2] = normal[i] / normal_length;
        }
        return m;
    }

    void MultiplyMatrixVector3x3(coord &y, const matrix &m, const coord &x) {
        for (int i = 0; i < 3; ++i) {
            y[i] = 0;
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                y[i] += m[i][j] * x[j];
            }
        }
    }

    void MultiplyMatrixMatrix3x3(matrix& res, const matrix& m1, const matrix& m2) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                res[i][j] = 0;
            }
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    res[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
    }

    void MultiplyVectorNumber(coord &y, double num) {
        for (int i = 0; i < 3; ++i) {
            y[i] *= num;
        }
    }

    void PrintMatrix3x3(const matrix& a, std::ostream& out) {
        out << a[0][0] << ' ' << a[0][1] << ' ' << a[0][2] << '\n';
        out << a[1][0] << ' ' << a[1][1] << ' ' << a[1][2] << '\n';
        out << a[2][0] << ' ' << a[2][1] << ' ' << a[2][2] << '\n';
    }
}