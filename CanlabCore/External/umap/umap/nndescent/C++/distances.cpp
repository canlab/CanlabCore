//
// Created by Stephen Meehan on 11/12/20.
//
/*
AUTHORS
    Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
    Stephen Meehan <swmeehan@stanford.edu>
    Connor Meehan <connor.gw.meehan@gmail.com>

Provided by suh ( Stanford University's Herzenberg Lab)
License: BSD 3 clause
*/

#include "distances.h"
#include <cmath>

using namespace suh;
using std::abs; using std::pow;


class Euclidean : public DistanceMetric {
public:
    double compute(const double *x, const double *y, const size_t size) override {
        double result = 0.0;
        for (const double *const end = x + size; x < end; x++, y++)
            result += ((*x - *y) * (*x - *y));
        return sqrt(result);
    }

    double finalize_sortable(const double sortable_distance) override {
        return sqrt(sortable_distance);
    }

    double compute_sortable(const double *x, const double *y, const int size) override {
        double result = 0.0;
        for (const double *const end = x + size; x < end; x++, y++)
            result += ((*x - *y) * (*x - *y));
        return result;
    }

    Euclidean() {
        needs_finalizing_ = true;
    }
};

class SEuclidean : public DistanceMetric {
public:
    suh::RangeSquaredMatrix<double> rangeSquared;

public:


    double compute(const double *x, const double *y, const size_t size) override {
        double result = 0.0;

        for (size_t j = 0; j < size; j++) {
            result += (x[j] - y[j]) * (x[j] - y[j]) / rangeSquared.matrix_[0][j];
        }

        return sqrt(result);
    }

    double finalize_sortable(const double sortable_distance) override {
        return sqrt(sortable_distance);
    }

    double compute_sortable(const double *x, const double *y, const int size) override {
        double result = 0.0;

        for (size_t j = 0; j < size; j++) {
            result += (x[j] - y[j]) * (x[j] - y[j]) / rangeSquared.matrix_[0][j];
        }

        return result;
    }


    SEuclidean( double **prior_matrix, const size_t columns)
            : rangeSquared(prior_matrix, columns){

    }

    SEuclidean(const double **data, const size_t rows, const size_t columns)
            : rangeSquared(data, rows, columns) {
        needs_finalizing_ = true;
    }
};

class CityBlock : public DistanceMetric {
public:
    double compute(const double *x, const double *y, const size_t size) override {
        double result = 0.0;
        for (const double *const end = x + size; x < end; x++, y++)
            result += abs(*x - *y);
        return result;
    }

    CityBlock() {
    }
};

class Chebychev : public DistanceMetric {
public:
    double compute(const double *x, const double *y, const size_t size) override {
        double result = 0.0;
        for (const double *const end = x + size; x < end; x++, y++)
            result = (result < abs(*x - *y)) ? abs(*x - *y) : result;
        return result;
    }

    Chebychev() {
    }
};

class Minkowski : public DistanceMetric {
public:
    double compute(const double *x, const double *y, const size_t size) override {
        double result = 0.0;
        for (const double *const end = x + size; x < end; x++, y++)
            result += pow(abs(*x - *y), p);
        return pow(result, 1. / p);
    }

    double finalize_sortable(const double sortable_distance) override {
        return pow(sortable_distance, 1. / p);
    }

    double compute_sortable(const double *x, const double *y, const int size) override {
        double result = 0.0;
        for (const double *const end = x + size; x < end; x++, y++)
            result += pow(abs(*x - *y), p);
        return result;
    }

    Minkowski() {
        needs_finalizing_ = true;
    }

    Minkowski(double exp) {
        needs_finalizing_ = true;
        p = exp;
    }

protected:
    double p = 2; // same default used by MatLab
};

class Mahalanobis : public DistanceMetric {
    suh::InverseCovarianceMatrix<double> invCov;

public:
    double compute(const double *x, const double *y, const size_t size) override {
        double result = 0.0;
        for (size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < size; j++) {
                result += invCov.matrix_[i][j] * (x[i] - y[i]) * (x[j] - y[j]);
            }
        }
        return sqrt(result);
    }

    double finalize_sortable(const double sortable_distance) override {
        return sqrt(sortable_distance);
    }

    double compute_sortable(const double *x, const double *y, const int size) override {
        double result = 0.0;
        for (size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < size; j++) {
                result += invCov.matrix_[i][j] * (x[i] - y[i]) * (x[j] - y[j]);
            }
        }
        return result;
    }


    Mahalanobis( double **const prior_matrix, const size_t columns)
            : invCov(prior_matrix, columns) {
    }

    Mahalanobis(const double **data, const size_t rows, const size_t columns) :
            invCov(data, rows, columns) {
        needs_finalizing_ = true;
    }
};

class Hamming : public DistanceMetric {
public:
    double compute(const double *x, const double *y, const size_t size) override {
        double result = 0.0;
        for (const double *const end = x + size; x < end; x++, y++)
            result += (*x == *y) ? 0. : 1.;
        return result / size;
    }

    Hamming() {
    }
};

class Jaccard : public DistanceMetric {
public:
    double compute(const double *x, const double *y, const size_t size) override {
        double difference = 0.0;
        double support = 0.0;

        for (const double *const end = x + size; x < end; x++, y++)
            if ((*x != 0.) || (*y != 0.)) {
                support++;
                difference += (*x == *y) ? 0. : 1.;
            }
        if (support == 0.)
            return 0.;
        else
            return difference / support;
    }

    Jaccard() {
    }
};

class Cosine : public DistanceMetric {
public:
    double compute(const double *x, const double *y, const size_t size) override {
        double xNormSquared = 0.0;
        double yNormSquared = 0.0;
        double dotProduct = 0.0;

        for (const double *const end = x + size; x < end; x++, y++) {
            xNormSquared += (*x) * (*x);
            yNormSquared += (*y) * (*y);
            dotProduct += (*x) * (*y);
        }

        return 1. - dotProduct / sqrt(xNormSquared * yNormSquared);
    }

    double finalize_sortable(const double sortable_distance) override {
        if (sortable_distance < 0)
            return 1. - sqrt(-1. * sortable_distance);
        else
            return 1. + sqrt(sortable_distance);
    }

    double compute_sortable(const double *x, const double *y, const int size) override {
        double xNormSquared = 0.0;
        double yNormSquared = 0.0;
        double negDotProduct = 0.0;

        for (const double *const end = x + size; x < end; x++, y++) {
            xNormSquared += (*x) * (*x);
            yNormSquared += (*y) * (*y);
            negDotProduct -= (*x) * (*y);
        }

        if (negDotProduct < 0)
            return -1. * negDotProduct * negDotProduct / xNormSquared / yNormSquared;
        else
            return negDotProduct * negDotProduct / xNormSquared / yNormSquared;

    }

    Cosine() {
        needs_finalizing_ = true;
    }
};

class Correlation : public DistanceMetric {
public:
    double compute(const double *x, const double *y, const size_t size) override {
        double xNormSquared = 0.0;
        double xSum = 0.0;
        double yNormSquared = 0.0;
        double ySum = 0.0;
        double dotProduct = 0.0;

        for (const double *const end = x + size; x < end; x++, y++) {
            xNormSquared += (*x) * (*x);
            xSum += *x;
            yNormSquared += (*y) * (*y);
            ySum += *y;
            dotProduct += (*x) * (*y);
        }

        return 1. - (size * dotProduct - xSum * ySum) /
                    sqrt((size * xNormSquared - xSum * xSum) * (size * yNormSquared - ySum * ySum));
    }

    double finalize_sortable(const double sortable_distance) override {
        if (sortable_distance < 0)
            return 1. - sqrt(-1. * sortable_distance);
        else
            return 1. + sqrt(sortable_distance);
    }

    double compute_sortable(const double *x, const double *y, const int size) override {
        double xNormSquared = 0.0;
        double xSum = 0.0;
        double yNormSquared = 0.0;
        double ySum = 0.0;
        double dotProduct = 0.0;

        for (const double *const end = x + size; x < end; x++, y++) {
            xNormSquared += (*x) * (*x);
            xSum += *x;
            yNormSquared += (*y) * (*y);
            ySum += *y;
            dotProduct += (*x) * (*y);
        }

        double numerator = xSum * ySum - size * dotProduct;

        if (numerator < 0)
            return -1. * numerator * numerator / (size * xNormSquared - xSum * xSum) /
                   (size * yNormSquared - ySum * ySum);
        else
            return numerator * numerator / (size * xNormSquared - xSum * xSum) / (size * yNormSquared - ySum * ySum);

    }

    Correlation() {
        needs_finalizing_ = true;
    }
};


DistancePtr DistanceMetric::euclidean() {
    static DistancePtr singleton(new Euclidean());
    return singleton;
}

#define IS(r) suh::strcmpi(name,r)==0

bool DistanceMetric::IsSupported(const char *name) {
    bool supported = true;
    if (IS("euclidean")) {
    } else if (IS("cityblock")) {
    } else if (IS("cosine")) {
    } else if (IS("mahalanobis")) {
    } else if (IS("seuclidean")) {
    } else if (IS("minkowski")) {
    } else if (IS("chebychev")) {
    } else if (IS("correlation")) {
    } else if (IS("hamming")) {
    } else if (IS("jaccard")) {
    } else if (IS("spearman")) {
    } else {
        supported = false;
    }
    return supported;
}

DistancePtr DistanceMetric::find(const std::string &name,
                                 const void *dist_args, const double **self,
                                 const size_t self_rows, const size_t self_columns,
                                 const double **other, const size_t other_rows,
                                 const size_t other_columns) {
    bool unrecognized = false;
    DistanceMetric *result = nullptr;
    //bool ok = IsSupported(name.c_str());
    if (IS("euclidean")) {
        result = new Euclidean();
    } else if (IS("cityblock")) {
        result = new CityBlock();
    } else if (IS("cosine")) {
        result = new Cosine();
    } else if (IS("mahalanobis")) {
        if (dist_args == nullptr)
            result = new Mahalanobis(other, other_rows, other_columns);
        else
            result = new Mahalanobis((double **)dist_args, self_columns);
    } else if (IS("seuclidean")) {
        if (dist_args==nullptr)
            result = new SEuclidean(self, self_rows, self_columns);
        else
            result = new SEuclidean((double **) dist_args, self_columns);
    } else if (IS("minkowski")) {
        if (dist_args == nullptr)
            result = new Minkowski();
        else
            result=new Minkowski(*((double *)dist_args));
    } else if (IS("chebychev")) {
        result = new Chebychev();
    } else if (IS("correlation")) {
        result = new Correlation();
    } else if (IS("hamming")) {
        result = new Hamming();
    } else if (IS("jaccard")) {
        result = new Jaccard();
    } else if (IS("spearman")) {
        // TODO: do equivalent work to "cityblock"
    } else {
        unrecognized = true;
        std::cerr << "UNRECOGNIZED metric name \"" << name << "\"!!!" << std::endl;
    }
    if (result == nullptr) {
        if (!unrecognized)
            std::cerr << "SIGH ...implementation coming soon for \"" << name << "\" ... " << std::endl;
        result = new Euclidean();
        std::cerr << "defaulting to \"euclidean\"" << std::endl;
    }
    result->name = name;
    return std::shared_ptr<DistanceMetric>(result);
}


DistanceMetric::~DistanceMetric() {
    if (suh::debug_ctor_dtor) {
        std::cout << "Destructing distance for " << name << std::endl;
    }
}

