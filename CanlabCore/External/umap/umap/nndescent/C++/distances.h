//
// Created by Jonathan Ebrahimian on 7/17/20.
//
/*

AUTHORS
    Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
    Connor Meehan <connor.gw.meehan@gmail.com>
    Stephen Meehan <swmeehan@stanford.edu>

Provided by suh ( Stanford University's Herzenberg Lab)
License: BSD 3 clause
*/

// Expect to evolve this file when we start adding other distance metrics
// Also plan to deprecate preprocessor macros (e.g. euclideanMacro)
#ifndef SUH_DISTANCES_H
#define SUH_DISTANCES_H

#include "suh_math.h"

namespace suh {
    inline double euclidean( const double * x,  const double *  y, const size_t size){
        double result = 0.0;
        // accelerate slightly with pointer arithmetic + indirection
        // instead of integer arithmetic + indexing+indirection
        for (const double * const end=x+size;x<end;x++, y++)
            result +=((*x - *y)*(*x - *y));
        return sqrt(result);
    }

    // DistanceMetric is an abstract class for computing one of the following types of distances
/*
     *  'euclidean'         Euclidean distance (default).
     *  'cityblock'         City block distance.
     *  'cosine'            One minus the cosine of the included angle between points (treated as vectors).
     *  'mahalanobis'       Mahalanobis distance using the sample covariance of X, C = cov(X,'omitrows').
     *                      Use sub class data to specify another value for C, where the matrix C is
     *                      symmetric and positive definite.
     *  'squaredeuclidean'  Squared Euclidean distance. (This option is provided for efficiency only.
     *                      It does not satisfy the triangle inequality.)
     *  'seuclidean'        Standardized Euclidean distance. Each coordinate difference between observations
     *                      is scaled by dividing by the corresponding element of the standard deviation,
     *                      S = std(X,'omitnan'). Use sub class data to specify another value for S.
     *  'minkowski'         Minkowski distance. The default exponent is 2. Use sub class data to specify a
     *                      different exponent P, where P is a positive scalar value of the exponent.
     *  'chebychev'         Chebychev distance (maximum coordinate difference).
     *  'correlation'       One minus the sample correlation between points (treated as sequences of values).
     *  'hamming'           Hamming distance, which is the percentage of coordinates that differ.
     *  'jaccard'           One minus the Jaccard coefficient, which is the percentage of nonzero coordinates that differ.
     *  'spearman'          One minus the sample Spearman's rank correlation between observations
     *                      (treated as sequences of values).
*/

    class DistanceMetric {
        std::string name;

    public:
        virtual double compute(const double * x,  const double *  y, const size_t  size)=0;
        virtual double compute_sortable(const double * x,  const double *  y, const int  size) {
            // if no sortable version exists then just use full computation
            return compute(x,y,size);
        }

        // for efficiency this computes an intermediate form of the distance for sorting purposes with KnnGraph
        virtual double finalize_sortable(const double sortable_distance) {
            // if no sortable version exists then just use full computation
            return sortable_distance;
        }

        inline static std::shared_ptr<DistanceMetric> find(std::string name,
                               const void * dist_args, const suh::MatrixPtr self) {
            return find(name.c_str(), dist_args, (const double **) self->matrix_, self->rows_, self->columns_);
        }

        inline static std::shared_ptr<DistanceMetric> find(std::string name,
                           const void *dist_args, const suh::MatrixPtr self, const suh::MatrixPtr other){
            return find(name.c_str(), dist_args, (const double **)self->matrix_, self->rows_, self->columns_,
                        (const double **)other->matrix_, other->rows_, other->columns_);
        }
        static std::shared_ptr<DistanceMetric>find(const std::string &name, const void *dist_args,
               const double** self, const size_t self_rows, const size_t self_columns,
               const double ** other=nullptr, const size_t other_rows=0, const size_t other_columns=0);


        virtual ~DistanceMetric();

        static std::shared_ptr<DistanceMetric> euclidean();

    protected:
        // set by sub class if sortable distance is worth computing separately
        bool needs_finalizing_=false;
    public:
        bool needs_finalizing(){
            return needs_finalizing_;
        }
        void set_finalizing_complete(){
            needs_finalizing_=false;
        }

        static bool IsSupported(const char *name);
    };

    using DistancePtr = std::shared_ptr<DistanceMetric>;
}

#endif //SUH_DISTANCES_H
