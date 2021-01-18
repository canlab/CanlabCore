//
// Created by Jonathan Ebrahimian on 7/16/20.
//
/*
AUTHORS
 Original Python:
    Leland McInnes <leland.mcinnes@gmail.com>
 C++ translation & optimization:
    Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
    Connor Meehan <connor.gw.meehan@gmail.com>
    Stephen Meehan <swmeehan@stanford.edu>

ALGORITHMS
   1 Dong, Wei, Charikar, Moses, Li, Kai; 
     Efficient K-Nearest Neighbor Graph Construction for Generic Similarity Measures; 
     https://www.cs.princeton.edu/cass/papers/www11.pdf 
   2 Dasgupta, Sanjoy, and Freund, Yoav; 
     Random projection trees and low dimensional manifolds;
     https://cseweb.ucsd.edu/~dasgupta/papers/rptree-stoc.pdf 

Provided by suh ( Stanford University's Herzenberg Lab)
License: BSD 3 clause
*/

#include "KnnGraph.h"
#include "suh.h"
#include <cmath>
#include <iostream>
#include <limits>

#include <iomanip>
#include <string>

using namespace suh;

KnnGraph::KnnGraph(const int rows, const int cols, DistancePtr dist, const bool fillNow)
        : nn_graph_(nullptr), rows_(rows), columns_(cols), dist(dist) {
    if (fillNow) {
        double defaultValues[] = {-1, INFINITY, 0};
        nn_graph_ = suh::new_matrices<double>(3, rows, cols, defaultValues);
    } else {
        nn_graph_ = suh::new_matrices<double>(3, rows, cols);
    }
    if (!dist) {
        this->dist = DistanceMetric::euclidean();
    }
}

KnnGraph::~KnnGraph() {
    suh::delete_matrices(nn_graph_, 3);
    if (suh::debug_ctor_dtor) {
        std::cout << "Destructing " << (std::string) *this;
        if (id >= 0)
            std::cout << " id=" << id;
        if (nn_graph_ == nullptr) {
            std::cout << " (previously deallocated)";
        }
        std::cout << std::endl;
    }
}

KnnGraph::KnnGraph(const KnnGraph &other)
        : nn_graph_(other.nn_graph_), rows_(other.rows_), columns_(other.columns_) {
    if (suh::debug_ctor_dtor)
        std::cout << "Copy constructor from " << other.id << " does work of allocating and copying" << std::endl;
    nn_graph_ = new_matrices<double>((const double ***) other.nn_graph_, 3, rows_, columns_);
}

KnnGraph::KnnGraph(double ***prior, const int rows, const int columns, DistancePtr dist, const bool copy)
        : nn_graph_(prior), rows_(rows), columns_(columns), dist(dist) {
    if (copy){
        nn_graph_ = suh::new_matrices<double>(3, rows, columns);
        for (int i=0;i<3;i++){
            suh::copy(nn_graph_[i], (const double **)prior[i], rows, columns);
        }
    }
}


KnnGraph::KnnGraph(KnnGraph &&other)
        : nn_graph_(other.nn_graph_), rows_(other.rows_), columns_(other.columns_), dist(other.dist) {
    if (suh::debug_ctor_dtor)
        std::cout << "Move  constructor from " << other.id << " just transfers rvalue pointers and sets to nullptr"
                  << std::endl;
    other.nn_graph_ = nullptr;
}


KnnGraph &KnnGraph::operator=(const KnnGraph &other) {
    this->~KnnGraph();
    new(this)KnnGraph(other);
    return *this;
}

KnnGraph &KnnGraph::operator=(KnnGraph &&other) {
    this->~KnnGraph();
    new(this)KnnGraph(std::forward<KnnGraph>(other));
    return *this;
}

MatrixPtr KnnGraph::indices() const {
    Matrix<double> inds(nn_graph_[0][0], rows_, columns_, Transfer::kCopy);
    return inds.shared_ptr();
}

void KnnGraph::copy_indices_column_wise(double *to) const {
    Matrix<double> t(nn_graph_[0], rows_, columns_);
    t.copy_column_wise(to);
    // limit cost of copying to table
    t.block_matrix_deallocation();
}

void KnnGraph::copy_indices_MatLab(double *to) const {
    Matrix<double> t(nn_graph_[0], rows_, columns_);
    t.copy_column_wise(to, 1);
    // limit cost of copying to table
    t.block_matrix_deallocation();
}

void KnnGraph::copy_distances_column_wise(double *to) const {
    Matrix<double> t(nn_graph_[1], rows_, columns_);
    t.copy_column_wise(to);
    // limit cost of copying to table
    t.block_matrix_deallocation();
}

MatrixPtr KnnGraph::distances() const {
    if (dist->needs_finalizing()){
        for (double *start=nn_graph_[1][0], *end=start+(rows_*columns_);
                     start<end; start++)
            *start=dist->finalize_sortable(*start);
        dist->set_finalizing_complete();
    }
    Matrix<double> dists(nn_graph_[1][0], rows_, columns_, Transfer::kCopy);
    return dists.shared_ptr();
}

KnnGraph KnnGraph::debug_move_copy(const KnnGraph &in) {
    //in.id=99;
    return in;
}

KnnGraph KnnGraph::debug_move_copy(double ***ptr, const int rows, const int columns) {
    KnnGraph graphCopy1(new_matrices((const double ***) ptr, 3, rows, columns),
                        rows, columns, DistanceMetric::euclidean(), false);
    graphCopy1.id = 22; // identify for destructing
    KnnGraph graphCopy2(graphCopy1);
    graphCopy2.id = 33;
    KnnGraph graphMove(std::move(graphCopy1));
    graphMove.id = 44; // identify for destructing
    graphCopy1 = graphCopy2;
    graphCopy1.id = 55;
    graphMove = std::move(graphCopy2);
    graphMove.id = 66;
    double ***ptr2 = new_matrices((const double ***) ptr, 3, rows / 2, columns / 4);
    //KnnGraph graphMove2=debug_move_copy(KnnGraph(ptr2, rows/2, columns/3));
    KnnGraph graphMove3 = debug_move_copy(std::move(KnnGraph(ptr2, rows / 2, columns / 3, DistanceMetric::euclidean(), false)));
    return graphMove;
}

// (Connor) Seems to be necessary for Windows computers
void KnnGraph::init(double ***const nn_graph, const int rows, const int columns) {
    double initialValue;
    for (int i = 0; i < 3; i++) {
        if (i == 0) {
            initialValue = -1;
        } else if (i == 1) {
            initialValue = INFINITY;
        } else {
            initialValue = 0;
        }
        std::fill_n(nn_graph[i][0], rows * columns, initialValue);
    }
}

int KnnGraph::push(double ***const nn_graph, const int row,
                   const double dist, // original Python named this weight even though the heap was for nearest neighbor DISTANCE
                   const int index, const int flag, const int columns) {
    double *dists = nn_graph[1][row];
    if (dist >= nn_graph[1][row][0]) {
        return 0;
    }
    const double dIndex = index;// avoid constant cost of coverting int to double
    double *indices = nn_graph[0][row], *is_new = nn_graph[2][row];
    for (double *start = indices, *end = indices + columns; start != end; start++) {
        if (dIndex == *start) {
            return 0;
        }
    }
// insert val at position zero
    *dists = dist;
    *indices = dIndex;
    *is_new = flag;
    int i = 0;
    int i_swap;
    while (true) {
        int ic1 = 2 * i + 1;
        int ic2 = ic1 + 1;

        if (ic1 >= columns) {
            break;
        } else if (ic2 >= columns) {
            if (dists[ic1] > dist) {
                i_swap = ic1;
            } else {
                break;
            }
        } else if (dists[ic1] >= dists[ic2]) {
            if (dist < dists[ic1]) {
                i_swap = ic1;
            } else {
                break;
            }
        } else {
            if (dist < dists[ic2]) {
                i_swap = ic2;
            } else {
                break;
            }
        }

        dists[i] = dists[i_swap];
        indices[i] = indices[i_swap];
        is_new[i] = is_new[i_swap];

        i = i_swap;
    }
    dists[i] = dist;
    indices[i] = dIndex;
    is_new[i] = flag;
    return 1;
}

int
KnnGraph::unchecked_push(double ***const nn_graph, const int row, const double dist, const int index, const int flag,
                         const int columns) {
    double *dists = nn_graph[1][row];
    if (dist >= nn_graph[1][row][0]) {
        return 0;
    }
    const double dIndex = index;// save on constant conversion of int to double
    double *indices = nn_graph[0][row],
            *is_new = nn_graph[2][row];

// insert val at position zero
    *dists = dist;
    *indices = dIndex;
    *is_new = flag;

    int i = 0;
    int i_swap;
    while (true) {
        int ic1 = 2 * i + 1;
        int ic2 = ic1 + 1;

        if (ic1 >= columns) {
            break;
        } else if (ic2 >= columns) {
            if (dists[ic1] > dist) {
                i_swap = ic1;
            } else {
                break;
            }
        } else if (dists[ic1] >= dists[ic2]) {
            if (dist < dists[ic1]) {
                i_swap = ic1;
            } else {
                break;
            }
        } else {
            if (dist < dists[ic2]) {
                i_swap = ic2;
            } else {
                break;
            }
        }

        dists[i] = dists[i_swap];
        indices[i] = indices[i_swap];
        is_new[i] = is_new[i_swap];

        i = i_swap;
    }
    dists[i] = dist;
    indices[i] = dIndex;
    is_new[i] = flag;
    return 1;
}

void KnnGraph::rejection_sample(const int columns, const int rows, long *const rng_state, int *const indices) {
    int j = 0;
    for (int i = 0; i < columns; i++) {
        bool flagFound = false;
        while (true) {
            j = tau_rand_int(rng_state, rows);
            for (int *start = indices, *end = indices + i; start != end; start++) {
                if (j == *start) {
                    flagFound = true;
                    break;
                }
            }
            if (!flagFound) {
                break;
            } else {
                flagFound = false;
            }
        }
        indices[i] = j;
    }
}

void KnnGraph::new_build_candidates(
        double ***const nn_graph, const int rows, const int columns, const int max_candidates,
        long *const rng_state, double ***const new_candidate_neighbors, double ***const old_candidate_neighbors,
        const double rho) {
    KnnGraph::init(new_candidate_neighbors, rows, max_candidates);
    KnnGraph::init(old_candidate_neighbors, rows, max_candidates);

    int idx;
    int isn;
    double d;
    int c;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (nn_graph[0][i][j] < 0) {
                continue;
            }
            idx = nn_graph[0][i][j];
            isn = nn_graph[2][i][j];
            d = TauRand(rng_state);

            if (TauRand(rng_state) < rho) {
                c = 0;
                if (isn != 0) {
                    c += push(new_candidate_neighbors, i, d, idx, isn, max_candidates);
                    c += push(new_candidate_neighbors, idx, d, i, isn, max_candidates);
                    if (c > 0) {
                        nn_graph[2][i][j] = 0;
                    }
                } else {
                    push(old_candidate_neighbors, i, d, idx, isn, max_candidates);
                    push(old_candidate_neighbors, idx, d, i, isn, max_candidates);
                }
            }
        }
    }
    return;
}


/*TODO: look at sift down implementation. Could it be optimized?
 *
 *
 */
void KnnGraph::sort(double ***const nn_graph, const int rows, const int columns) {

    double temp;
    double temp2;

    double *tempVec = new double[columns];
    double *tempVec2 = new double[columns];

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns - 1; j++) {
            temp = nn_graph[0][i][0];
            temp2 = nn_graph[0][i][columns - j - 1];
            nn_graph[0][i][0] = temp2;
            nn_graph[0][i][columns - j - 1] = temp;

            temp = nn_graph[1][i][0];
            temp2 = nn_graph[1][i][columns - j - 1];
            nn_graph[1][i][0] = temp2;
            nn_graph[1][i][columns - j - 1] = temp;
            int cutoff = columns - j - 1;
            slice(nn_graph[1][i], cutoff, tempVec);
            slice(nn_graph[0][i], cutoff, tempVec2);
            sift_down(tempVec, tempVec2, cutoff);
            for (int x = 0; x < cutoff; x++) {
                nn_graph[1][i][x] = tempVec[x];
                nn_graph[0][i][x] = tempVec2[x];
            }
        }
    }
    delete[] tempVec;
    delete[] tempVec2;
}


int KnnGraph::smallest_flagged(double ***const nn_graph, const int row, const int columns) {
    double *indices = nn_graph[0][row], *dists = nn_graph[1][row], *is_new = nn_graph[2][row];

    double min_dist = std::numeric_limits<double>::infinity();
    int result_index = -1;

    for (int i = 0; i < columns; i++) {
        if (is_new[i] == 1 && dists[i] < min_dist) {
            min_dist = dists[i];
            result_index = i;
        }
    }
    if (result_index >= 0) {
        is_new[result_index] = 0.0;
        return int(indices[result_index]);
    } else {
        return -1;
    }
}
