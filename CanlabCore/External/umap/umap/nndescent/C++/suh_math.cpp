//
// Created by Stephen Meehan on 12/8/20.
//
#include "suh_math.h"

namespace suh {



//Find the cofactor matrix of a square matrix
    void Cofactor( double **a, double **b, const size_t n) {
        int i, j, ii, jj, i1, j1;
        double det;
        double **c;

        c = static_cast<double **>(malloc((n - 1) * sizeof(double *)));
        for (i = 0; i < n - 1; i++)
            c[i] = static_cast<double *>(malloc((n - 1) * sizeof(double)));

        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {

                /* Form the adjoint a_ij */
                i1 = 0;
                for (ii = 0; ii < n; ii++) {
                    if (ii == i)
                        continue;
                    j1 = 0;
                    for (jj = 0; jj < n; jj++) {
                        if (jj == j)
                            continue;
                        c[i1][j1] = a[ii][jj];
                        j1++;
                    }
                    i1++;
                }

                /* Calculate the determinate */
                det = determinant(c, n - 1);

                /* Fill in the elements of the cofactor */
                b[i][j] = pow(-1.0, i + j) * det;
            }
        }
        for (i = 0; i < n - 1; i++)
            free(c[i]);
        free(c);
    }


//   Transpose of a square matrix, do it in place
    void TransposeScale(double **a, const int n, const double scalar) {
        int i, j;
        double tmp;

        for (i = 1; i < n; i++) {
            for (j = 0; j < i; j++) {
                tmp = a[i][j];
                a[i][j] = a[j][i] * scalar;
                a[j][i] = tmp * scalar;
            }
        }

        for (i = 0; i < n; i++) {
            a[i][i] = a[i][i] * scalar;
        }
    }

// Function to calculate and store inverse, returns false if
// matrix is singular
    bool inverse(double **a, double **inv, const size_t n) {
        std::cout << "1st determinant computing now"<<std::endl;
        double det=determinant(a,n);
        //double dbg2  = Determinant(a, n);
        std::cout << "doing CoFactor now"<<std::endl;
        if (det == 0) {
            std::cout << "Singular matrix, can't find its inverse" << std::endl;
            return false;
        }

        Cofactor(a, inv, n);
        std::cout << "CoFactor done"<<std::endl;
        TransposeScale(inv, n, 1. / det);

        return true;
    }

    double Determinant(double **a,size_t n)
    {
        int i,j,j1,j2;
        double det = 0;
        double **matrix_ = NULL;

        if (n < 1) { /* Error */

        } else if (n == 1) { /* Shouldn't get used */
            det = a[0][0];
        } else if (n == 2) {
            det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
        } else {
            det = 0;
            for (j1=0;j1<n;j1++) {
                matrix_ = static_cast<double **>(malloc((n - 1) * sizeof(double *)));
                for (i=0;i<n-1;i++)
                    matrix_[i] = static_cast<double *>(malloc((n - 1) * sizeof(double)));
                for (i=1;i<n;i++) {
                    j2 = 0;
                    for (j=0;j<n;j++) {
                        if (j == j1)
                            continue;
                        matrix_[i - 1][j2] = a[i][j];
                        j2++;
                    }
                }
                det += pow(-1.0,j1) * a[0][j1] * Determinant(matrix_, n - 1);
                for (i=0;i<n-1;i++)
                    free(matrix_[i]);
                free(matrix_);
            }
        }
        return(det);
    }

    class Determinator {
    public:

        const size_t  N;
         double *pw;
        Matrix<double> **pool;
        Determinator(const size_t N) : N(N){
            pool=new Matrix<double>*[N];
            pw=new double[N];
            for (int i=0;i<N;i++){
                pw[i]=pow(-1.0,i);
                int n=N-i;
                if (n>2){
                    pool[n-1]=new Matrix<double>(n-1,n-1);
                } else {
                    pool[n-1]= nullptr;
                }
            }
        }

        ~Determinator(){
            for (int i=0;i<N;i++){
                delete pool[i];
            }
            delete []pool;
        }
        double compute(const double **a,const size_t n);
    };

    double Determinator::compute(const double **a, const size_t n){
        int i,j,j1,j2;
        double det = 0;

        if (n < 1) { /* Error */
        } else if (n == 1) { /* Shouldn't get used */
            det = a[0][0];
        } else if (n == 2) {
            det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
        } else {
            Matrix<double> *matrix=pool[n-1];
            double **matrix_ = matrix->matrix_;

            det = 0;

            for (j1=0;j1<n;j1++) {
                //matrix->clear();
                for (i=1;i<n;i++) {
                    j2 = 0;
                    for (j=0;j<n;j++) {
                        if (j == j1)
                            continue;
                        matrix_[i - 1][j2] = a[i][j];
                        j2++;
                    }
                }
                det += pw[j1] * a[0][j1] * compute((const double **)matrix_, n - 1);
            }
        }
        return(det);
    }

    double determinant(const double **a, const size_t n){
        Determinator d=Determinator(n);
        const double result=d.compute(a,n);
        return result;
    }
}