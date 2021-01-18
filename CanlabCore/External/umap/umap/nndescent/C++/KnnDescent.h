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


#ifndef SUH_NNDESCENT_H
#define SUH_NNDESCENT_H
#include <vector>
#include <unordered_set>
// Not sure IF pair_hash is still needed
//#include "pair_hash.h"
#include "suh.h"
#include "KnnGraph.h"

namespace suh{


class KnnDescent {

public:
    static int numHeapPush;
    static int numDistances;
    static int numUncheckedHeapPush;
    static  KnnGraphPtr search(
            MatrixPtr self, const int n_neighbors, long *const rng_state,
            std::string dist_metric="euclidean", const int n_async_tasks=3,
             FncProgress verbose= nullptr, void *dist_args= nullptr, const int max_candidates= 50,
             const int n_iters= 10, const double delta= .001, const double rho= .5,
            const bool leaf= true, const bool low_memory= true,
            int n_trees= -1, const bool angular=false);
    static  KnnGraphPtr search(
            suh::MatrixPtr otherData, suh::MatrixPtr selfData,
            suh::MatrixIntPtr otherIndptr, suh::MatrixIntPtr otherIndices,
            const int n_neighbors,  const float transform_queue_size,
            long *const rng_state,  std::string dist_metric="euclidian",
            const int n_async_tasks=3, FncProgress verbose= nullptr, void *dist_args= nullptr,
            const long angular= false);

private:
    static  KnnGraphPtr nn_descent( // in Leland McInnes Python this was the main function for nn_descent on 1 dataset
            MatrixPtr self, const int n_neighbors, long *const rng_state, std::vector<std::vector<int> >&leaf_array,
            const int max_candidates, const int n_iters, const double delta, const double rho, const bool leaf,
            const bool low_memboery, FncProgress fnc_progress, const int n_async_tasks,
            std::string dist_metric, const void *dist_args);

    static void init_rp_tree(double ** const, double *** const, std::vector<std::vector<int> > &, const int, const int, DistancePtr dist);
    static void init_rp_tree_high_memory(double ** const , double *** const , std::vector<std::vector<int> > & ,
                                         const int , const int , bool ** const, DistancePtr dist);
    static bool nn_descent_internal_low_memory(
            double ***const nn_graph, const double ** const data, const int rows, const int columns,
            const int n_neighbors, long *rng_state, const int max_candidates, const int n_iters,
            const double delta, const double rho, FncProgress verbose, DistancePtr dist);
    static bool nn_descent_internal_low_memory_async(
            double ***const nn_graph, const double ** const data, const int rows, const int columns,
            const int n_neighbors, long *rng_state, const int max_candidates, const int n_iters,
            const double delta, const double rho, FncProgress verbose, DistancePtr dist, const int n_async_tasks);

    static bool nn_descent_internal_high_memory(double *** const, const double ** const, const int , const int , const int , long *const, bool ** const, const int , const int , const double, const double , FncProgress, DistancePtr);
    static void init_from_random(const int, const double **const, const int, const int, const double **const, const int, double *** const, long *const, DistancePtr dist);
    static void init_from_tree(const RpTree *const, const double** const, const int, const int, const int, const double** const, const int, double*** const, long* const, DistancePtr dist);

    static double*** initialise_search(const std::vector<RpTree *>&, const double ** const, const int, const int, const double **const, const int, const int, long *const, DistancePtr dist);
    static bool initialized_nnd_search(
            double ***const nn_graph, const int nn_graph_columns,
            double **const otherData, const int *const indptr, const int *const indices,
            double ** const selfData, const int self_rows, const int self_columns,
            DistancePtr dist, const int n_async_tasks=3, FncProgress verbose= nullptr);

    static bool initialized_nnd_search(
            double ***const nn_graph, const int nn_graph_columns,
            double **const other_data, const int *const indptr, const int *const indices,
            double **const self_data, const int self_rows, const int self_columns,
            DistancePtr dist, FncProgress verbose, const int start, const int end);

    class OtherSearchAsyncTask {
        friend KnnDescent;
        double *** const nn_graph;
        const int nn_graph_columns;
        double **const otherData;
        const int *const indptr;
        const int *const indices;
        double **const selfData;
        DistancePtr dist;
        const int selfRows, selfColumns, start, end;
        FncProgress  verbose;
        OtherSearchAsyncTask(double *** const nn_graph, const int nn_graph_columns,
                             double **const otherData, const int *const indptr, const int *const indices,
                             double **const selfData, const int selfRows, const int selfColumns,
                             DistancePtr dist, FncProgress verbose, const int start, const int end);
    public:
        void operator ()();
    };

    class SelfSearchAsyncTask {
        friend KnnDescent;
        double *** const nn_graph;
        const double ** const data;
        const int rows;
        const int columns;
        const int n_neighbors;
        long *const rng_state;
        const int max_candidates;
        const int n_iters;
        const double delta;
        const double rho;
        const int start;
        const int end;
        double ***old_candidates;
        double ***new_candidates;
        DistancePtr dist;
        SelfSearchAsyncTask(double ***const nnGraph, const double **const data, const int rows, const int columns,
                            const int n_neighbors, long *const rng_state, const int max_candidates, const int n_iters,
                            const double delta, const double rho, DistancePtr dist,
                            double ***old_candidates, double ***new_candidates, const int start, const int end);
    public:
        int operator()();
    };
};

}
#endif //SUH_NNDESCENT_H
