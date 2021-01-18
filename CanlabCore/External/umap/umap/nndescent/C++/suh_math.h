//
// Created by Stephen Meehan on 12/8/20.
//

#ifndef C___SUH_MATH_H
#define C___SUH_MATH_H

#include <cmath>
#include "suh.h"

namespace suh {

    //Find the cofactor matrix of a square matrix
    void Cofactor(double **a, double **b, const size_t n);

    //Transpose of a square matrix, do it in place
    void TransposeScale(double **a, const int n, const double scalar);

// Function to calculate and store inverse, returns false if
// matrix is singular
    bool inverse(double **a, double **inv, size_t n);

//  Recursive definition of determinant using expansion by minors.
    double Determinant(double **a, size_t n);

    double determinant(const double **a, const size_t n);

    inline double determinant(double **a, const size_t n) {
        return determinant((const double **) a, n);
    }

    template<class T>
    class InverseCovarianceMatrix : public Matrix<T> {
    public:
        InverseCovarianceMatrix(double **prior_matrix, const size_t columns) :
                Matrix<double>(prior_matrix, columns, columns) {
            this->block_matrix_deallocation();
        }

        InverseCovarianceMatrix(const T **data, const size_t rows, const size_t columns) : Matrix<double>(columns) {
            double *columnMeans = new double[columns];
            Matrix<double> cov(columns);

            for (size_t j = 0; j < columns; j++) {
                double sum = 0.0;
                for (size_t row = 0; row < rows; row++) {
                    sum += data[row][j];
                }
                columnMeans[j] = sum / rows;
            }

            for (size_t i = 0; i < columns; i++) {
                for (size_t j = i; j < columns; j++) {
                    double sum = 0.0;
                    for (size_t row = 0; row < rows; row++) {
                        sum += (data[row][i] - columnMeans[i]) * (data[row][j] - columnMeans[j]);
                    }
                    cov.matrix_[i][j] = sum / (rows - 1);
                }
            }
            delete[]columnMeans;

            for (size_t i = 0; i < columns; i++) {
                for (size_t j = 0; j < i; j++) {
                    cov.matrix_[i][j] = cov.matrix_[j][i];
                }
            }
            inverse(cov.matrix_, this->matrix_, columns);
        }
    };

    template<class T>
    class RangeSquaredMatrix : public Matrix<T> {
    public:

        RangeSquaredMatrix( double **prior_matrix, const size_t columns)
                : Matrix<double>(prior_matrix, 1, columns) {
            this->block_matrix_deallocation();
        }

        RangeSquaredMatrix(const T **data, const size_t rows, const size_t columns)
                : Matrix<double>(1, columns) {
            for (size_t j = 0; j < columns; j++) {
                double max = data[0][j];
                double min = data[0][j];
                for (size_t i = 1; i < rows; i++) {
                    max = (max < data[i][j]) ? data[i][j] : max;
                    min = (min > data[i][j]) ? data[i][j] : min;
                }
                if (max == min) {
                    std::cerr << "Range is 0 in one component!" << std::endl;
                    return;
                }
                this->matrix_[0][j] = (max - min) * (max - min);
            }
        }
    };

}
#endif //C___SUH_MATH_H
