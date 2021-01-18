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

/*
 * Use standard arrays instead of vectors
 * pass all vars by reference instead of by value
 */

#include "suh.h"
#include "KnnDescent.h"
#include "KnnGraph.h"
#include "RpTree.h"
#include "distances.h"
#include <chrono>
#include <iostream>
#include <unordered_set>
#include <set>
#include <cmath>
#include <iomanip>

#include <future>


using distFnc=double (*)(const double *, double *, const int size);

using namespace std::chrono;
using namespace suh;

int KnnDescent::numHeapPush;
int KnnDescent::numDistances;
int KnnDescent::numUncheckedHeapPush;

KnnGraphPtr KnnDescent::search(
        MatrixPtr self, const int n_neighbors, long *const rng_state,std::string dist_metric, const int n_async_tasks,
        FncProgress verbose, void *dist_args, const int max_candidates, const int n_iters, const double delta,
        const double rho,const bool leaf, const bool low_memory, int n_trees, const bool angular) {
    if (n_trees <= 0) {
        n_trees = 5 + int(round(pow((self->rows_), 0.5) / 20.0));
    }

    RpForestPtr rp_forest=RpForester::make_forest(
            self->matrix_, self->rows_, self->columns_, n_neighbors, n_trees,
            rng_state, angular);
    //rp_forest->test_serialization();

    std::vector<std::vector<int> > leaf_array = RpForester::leaf_array(rp_forest->get_grove());
    KnnGraphPtr knn_graph =KnnDescent::nn_descent(
            self, n_neighbors, rng_state, leaf_array, max_candidates, n_iters,
            delta, rho, leaf, low_memory, verbose, n_async_tasks, dist_metric, dist_args);
    if (knn_graph) {
        knn_graph->rp_forest_ = rp_forest;
    }
    return knn_graph;
}

KnnGraphPtr KnnDescent::nn_descent(
        MatrixPtr self, const int n_neighbors, long *const rng_state, std::vector<std::vector<int> > &leaf_array,
        const int max_candidates, const int n_iters, const double delta,
        const double rho, const bool rp_tree_init, const bool low_memory, FncProgress verbose,
        const int  n_async_tasks, std::string dist_metric, const void *dist_args) {
    DistancePtr dist=DistanceMetric::find(dist_metric, dist_args, self);

    bool** tried;
    const size_t columns=self->columns_;
    const size_t rows=self->rows_;
    double **data=self->matrix_;
    if(!low_memory){
        tried = new bool*[rows];
        for(int i = 0; i < rows; ++i){
            tried[i] = new bool[rows];
            for(int j = 0; j < rows; ++j){
                tried[i][j] = false;
            }
        }
    }
    KnnGraph *nnGraph=new KnnGraph(rows, n_neighbors, dist);
    double*** nn_graph = *nnGraph;
    int * indices = new int[n_neighbors];
    for(int i = 0; i < rows; i++){
        KnnGraph::rejection_sample(n_neighbors, rows, rng_state, indices);
        for(int j = 0; j < n_neighbors; j++){
            const double  d = dist->compute_sortable(data[i], data[indices[j]], columns);
            numDistances += 1;

            KnnGraph::push(nn_graph, i, d, indices[j], 1, n_neighbors);
            KnnGraph::push(nn_graph, indices[j], d, i, 1, n_neighbors);
            numHeapPush += 2;

            if(!low_memory){
                tried[i][indices[j]] = true;
                tried[indices[j]][i] = true;
            }
        }
    }
    delete [] indices;

    if(rp_tree_init && low_memory){
        init_rp_tree(data, nn_graph, leaf_array, n_neighbors, columns, dist);
    }else if(rp_tree_init && !low_memory){
        init_rp_tree_high_memory(data, nn_graph, leaf_array, n_neighbors, columns, tried, dist);
    }
    bool good=true;
    if(low_memory){
        if (n_async_tasks>1) {
            good=nn_descent_internal_low_memory_async(
                    nn_graph, (const double **) data, rows, columns, n_neighbors, rng_state, max_candidates, n_iters,
                    delta, rho, verbose, dist, n_async_tasks);
        } else{
            good=nn_descent_internal_low_memory(
                    nn_graph, (const double **) data, rows, columns, n_neighbors, rng_state, max_candidates, n_iters,
                    delta, rho, verbose, dist);
        }
    }else{
        good=nn_descent_internal_high_memory(
                nn_graph, (const double **)data, rows, columns, n_neighbors, rng_state, tried,
                max_candidates, n_iters, delta, rho,
                verbose, dist);
    }
    if (!good){
        return KnnGraphPtr();
    }
    KnnGraph::sort(nn_graph, rows, n_neighbors);
    return KnnGraphPtr(nnGraph);
}

void KnnDescent::init_rp_tree_high_memory(
        double ** const data, double *** const nn_graph, std::vector<std::vector<int> > &leaf_array,
        const int n_neighbors, const int columns, bool ** const tried, DistancePtr dist) {
    int p;
    int q;
    std::cout << "in high" << std::endl;
    for(int n = 0; n < leaf_array.size(); n++){
        for(int i = 0; i < leaf_array[n].size(); i++){
            p = leaf_array[n][i];
            if(p < 0){
                break;
            }
            for(int j = i + 1; j < leaf_array[n].size();j++){
                q = leaf_array[n][j];
                if(q < 0){
                    break;
                }
                if(tried[p][q]){
                    continue;
                }
                const double  d = dist->compute_sortable(data[p], data[q], columns);
                numDistances += 1;
                KnnGraph::push(nn_graph, p, d, q, 1, n_neighbors);
                numHeapPush += 1;
                tried[p][q] = true;
                if(p != q){
                    KnnGraph::push(nn_graph, p, d, q, 1, n_neighbors);
                    numHeapPush += 1;
                    tried[q][p] = true;
                }
            }
        }
    }
}

void KnnDescent::init_rp_tree(
        double ** const data, double *** const nn_graph, std::vector<std::vector<int> > &leaf_array,
        const int n_neighbors, const int columns, DistancePtr dist) {
    int p;
    int q;
    int push_in = 1;
    for(int n = 0; n < leaf_array.size(); n++){
        for(int i = 0; i < leaf_array[n].size(); i++){
            p = leaf_array[n][i];
            if(p < 0){
                break;
            }
            for(int j = i + 1; j < leaf_array[n].size(); j++){
                q = leaf_array[n][j];
                if(q < 0){
                    break;
                }
                const double  d = dist->compute_sortable(data[p], data[q], columns);
                numDistances += 1;
                KnnGraph::push(nn_graph, p, d, q, push_in, n_neighbors);
                numHeapPush += 1;
                if(p != q){
                    KnnGraph::push(nn_graph, q, d, p, push_in, n_neighbors);
                    numHeapPush += 1;
                }
            }
        }
    }
}


bool KnnDescent::nn_descent_internal_high_memory(
        double *** const nn_graph, const double **const data, const int rows, const int columns,
        const int n_neighbors, long * const rng_state, bool **const tried, const int max_candidates,
        const int n_iters, const double delta, const double rho, FncProgress verbose, DistancePtr dist) {
    int n_vertices = rows;
    KnnGraph oldGraph(rows, max_candidates, dist, false),newGraph(rows, max_candidates, dist, false);
    double*** new_candidate_neighbors = newGraph;
    double*** old_candidate_neighbors = oldGraph;
    int push_in = 1;

    int num = 0;
    for(int n = 0; n < n_iters; n++){
        if(verbose){
            if (!(verbose)(n, n_iters))
                return false;
        }
        KnnGraph::new_build_candidates(
                nn_graph, n_vertices, n_neighbors, max_candidates, rng_state,
                new_candidate_neighbors, old_candidate_neighbors, rho);

        int c = 0;
        int p;
        int q;
        long double d;

        for(int i = 0; i < n_vertices; i ++){
            for(int j = 0; j < max_candidates; j++){
                p = new_candidate_neighbors[0][i][j];
                if(p < 0){
                    num += 1;
                    continue;
                }

                for(int k = j; k < max_candidates; k++){
                    q = new_candidate_neighbors[0][i][k];

                    if (q < 0 || tried[p][q]){
                        continue;
                    }

                    const double  d = dist->compute_sortable(data[p], data[q], columns);
                    numDistances += 1;
                    c += KnnGraph::unchecked_push(nn_graph, p, d, q, push_in, n_neighbors);
                    numUncheckedHeapPush += 1;
                    tried[p][q] = true;
                    if(p != q){
                        c += KnnGraph::unchecked_push(nn_graph, q, d, p, push_in, n_neighbors);
                        numUncheckedHeapPush += 1;
                        tried[q][p] = true;
                    }
                }
                for(int k = 0;k < max_candidates;k++){
                    q = old_candidate_neighbors[0][i][k];
                    if(q < 0 || tried[p][q]){
                        continue;
                    }
                    const double  d = dist->compute_sortable(data[p], data[q], columns);
                    numDistances +=1 ;
                    c += KnnGraph::unchecked_push(nn_graph, p, d, q, push_in, n_neighbors);
                    numUncheckedHeapPush += 1;

                    tried[p][q] = true;
                    if(p != q){
                        c += KnnGraph::unchecked_push(nn_graph, q, d, p, push_in, n_neighbors);
                        numUncheckedHeapPush += 1;
                        tried[q][p] = true;
                    }
                }
            }
        }
        if(c <= delta * n_neighbors * rows){
            return true;
        }
    }
    return true;
}

bool KnnDescent::nn_descent_internal_low_memory(
        double *** const nn_graph, const double ** const data, const int rows, const int columns,
        const int n_neighbors, long * rng_state, const int max_candidates,
        const int n_iters, const double delta, const double rho, FncProgress verbose, DistancePtr dist){
    int push_in = 1;
    int num = 0;
    KnnGraph oldGraph(rows, max_candidates,  dist,false),newGraph(rows, max_candidates, dist, false);
    double*** new_candidate_neighbors = newGraph;
    double*** old_candidate_neighbors = oldGraph;

    for(int n = 0; n < n_iters; n++){
        if(verbose){
            if (!(verbose)(n, n_iters))
                return false;
        }
        KnnGraph::new_build_candidates(
                nn_graph, rows, n_neighbors, max_candidates, rng_state,
                new_candidate_neighbors, old_candidate_neighbors, rho);
        int c = 0;
        int p;
        int q;

        for(int i = 0; i < rows; i ++) {
            for (int j = 0; j < max_candidates; j++) {
                p = new_candidate_neighbors[0][i][j];
                if (p < 0) {
                    continue;
                }
                for(int k = j; k < max_candidates; k++){
                    q = new_candidate_neighbors[0][i][k];
                    if(q < 0){
                        num += 1;
                        continue;
                    }
                    const double  d = dist->compute_sortable(data[p], data[q], columns);
                    numDistances += 1;
                    c += KnnGraph::push(nn_graph, p, d, q, push_in, n_neighbors);
                    numHeapPush += 1;
                    if(p != q){
                        c += KnnGraph::push(nn_graph, q, d, p, push_in, n_neighbors);
                        numHeapPush += 1;
                    }
                }
                for(int k = 0; k < max_candidates; k++){
                    q = old_candidate_neighbors[0][i][k];
                    if(q < 0){
                        continue;
                    }
                    const double  d = dist->compute_sortable(data[p], data[q], columns);
                    numDistances += 1;
                    c += KnnGraph::push(nn_graph, p, d, q, push_in, n_neighbors);
                    numHeapPush += 1;
                    if(p != q){
                        c += KnnGraph::push(nn_graph, q, d, p, push_in, n_neighbors);
                        numHeapPush += 1;
                    }
                }
            }
        }
        if(c <= delta * n_neighbors * rows){
            return true;
        }
    }
    return true;
}

bool KnnDescent::nn_descent_internal_low_memory_async(
        double *** const nn_graph, const double ** const data, const int rows, const int columns,
        const int n_neighbors, long * rng_state, const int max_candidates,
        const int n_iters, const double delta, const double rho,
        FncProgress verbose, DistancePtr dist, const int n_async_tasks){
    int push_in = 1;
    KnnGraph oldGraph(rows, max_candidates, dist, false),newGraph(rows, max_candidates, dist, false);
    double*** new_candidate_neighbors = newGraph;
    double*** old_candidate_neighbors = oldGraph;
    for(int n = 0; n < n_iters; n++){
        if(verbose){
            if (!(verbose)(n, n_iters))
                return false;
        }
        KnnGraph::new_build_candidates(
                nn_graph, rows, n_neighbors, max_candidates, rng_state,
                new_candidate_neighbors, old_candidate_neighbors, rho);
        const int work= rows / n_async_tasks;
        std::vector<SelfSearchAsyncTask >tasks;
        std::vector<std::future<int>> futures;
        for (int i=0;i<n_async_tasks;i++) {
            const int start=i*work, end = i<n_async_tasks-1?start+work:rows;
            tasks.push_back(SelfSearchAsyncTask(
                    nn_graph, (const double **) data, rows, columns, n_neighbors,
                    rng_state, max_candidates, n_iters, delta, rho,  dist,
                    old_candidate_neighbors, new_candidate_neighbors, start, end));
        }
        for (SelfSearchAsyncTask &task:tasks) {
            auto f = std::async(task);
            futures.push_back(std::move(f));
        }
        int c = 0;
        for (auto &f:futures){
            c+=f.get();
        }
        if(c <= delta * n_neighbors * rows){
            return true;
        }
    }
    return true;
}

KnnGraphPtr KnnDescent::search(
        MatrixPtr otherData, MatrixPtr selfData,
        MatrixIntPtr otherIndptr, MatrixIntPtr otherIndices,
        const int n_neighbors, const float transform_queue_size,
        long *const rng_state, std::string dist_metric,
        const int n_async_tasks,FncProgress verbose,
        void *dist_args, const long angular){
    DistancePtr dist=DistanceMetric::find(dist_metric, dist_args, selfData, otherData);

    int n_trees = 5 + int(round(pow((otherData->rows_), 0.5) / 20.0));
    RpForestPtr rp_forest = RpForester::make_forest(otherData->matrix_, otherData->rows_, otherData->columns_, n_neighbors,
                                                    n_trees, rng_state, angular);
    //rp_forest->test_serialization();

    const int nn_graph_columns = (int) ((float)n_neighbors * transform_queue_size);
    double *** const nn_graph = KnnDescent::initialise_search(
            rp_forest->get_grove(), (const double **) otherData->matrix_, otherData->rows_, otherData->columns_,
            (const double **) selfData->matrix_, selfData->rows_, nn_graph_columns, rng_state, dist);
    KnnGraphPtr result;
    if (KnnDescent::initialized_nnd_search(
            nn_graph, nn_graph_columns,
            otherData->matrix_, otherIndptr->vector(), otherIndices->vector(),
            selfData->matrix_, selfData->rows_, selfData->columns_, dist, n_async_tasks, verbose)) {
        if (n_neighbors == transform_queue_size) {
            result = KnnGraphPtr(new KnnGraph(nn_graph, selfData->rows_, n_neighbors, dist, false));
        } else {
            result = KnnGraphPtr(new KnnGraph(nn_graph, selfData->rows_, n_neighbors, dist, true));
        }
        result->rp_forest_ = rp_forest;
    }
    return result;
}

double*** KnnDescent::initialise_search(
        const std::vector<RpTree *> & forest, const double ** const other_data, const int other_rows, const int columns,
        const double ** const self_data, const int self_rows, const int nn_graph_columns, long * const rng_state,
        DistancePtr dist) {
    KnnGraph *nnGraph=new KnnGraph(self_rows, nn_graph_columns);
    double***const nn_graph = *nnGraph;
    init_from_random(nn_graph_columns, other_data, other_rows, columns, self_data, self_rows, nn_graph, rng_state, dist);
    if(forest.size() != 0){
        for(int x = 0; x < forest.size(); x++){
            init_from_tree(forest[x], other_data, other_rows, columns, nn_graph_columns, self_data, self_rows, nn_graph,
                           rng_state, dist);
        }
    }
    return nn_graph;
}


void KnnDescent::init_from_random(
        const int nn_graph_columns, const double ** const other_data, const int other_rows, const int columns,
        const double ** const self_data, const int self_rows, double *** const nn_graph, long * const rng_state,
        DistancePtr dist){
    int * indices = new int[nn_graph_columns];
    for(int i = 0; i < self_rows; i++){
        KnnGraph::rejection_sample(nn_graph_columns, other_rows, rng_state, indices);
        for(int j = 0; j < nn_graph_columns; j++){
            if(indices[j] < 0){
                continue;
            }
            //std::cout << "figure out what size needs to go into euclidean calculation" << std::endl;
            const double  d = dist->compute_sortable(other_data[indices[j]], self_data[i], columns);
            numDistances += 1;
            KnnGraph::push(nn_graph, i, d, indices[j], 1, nn_graph_columns);
            numHeapPush += 1;
        }
        //suh::Debug(heap[0], i, 60);
    }
    delete [] indices;
    return;
}

void KnnDescent::init_from_tree(
        const RpTree *const tree, const double **const other_data,
        const int other_rows, const int columns,
        const int nn_graph_columns, const double **const self_data,
        const int self_rows, double ***const nn_graph, long *const rng_state, DistancePtr dist) {
    for (int i = 0; i < self_rows; i++ ){
        std::vector<int> indices = RpForester::search_flat_tree(self_data, i, columns, tree->hyperplanes, tree->offsets, tree->children, tree->indices, rng_state);//TODO implement
        for(int j = 0; j < indices.size(); j++){
            if(indices[j] < 0){
                continue;
            }
            const double  d = dist->compute_sortable(other_data[indices[j]], self_data[i], columns);
            numDistances += 1;
            KnnGraph::push(nn_graph, i, d, indices[j], 1, nn_graph_columns);
            numHeapPush += 1;
        }
    }
}


bool KnnDescent::initialized_nnd_search(
        double *** const nn_graph, const int nn_graph_columns,
        double **const otherData, const int *const indptr, const int *const indices,
        double **const selfData, const int self_rows, const int self_columns,
        DistancePtr dist, const int n_async_tasks, FncProgress verbose){
    int total_reports;
    if (n_async_tasks>1){
        total_reports=1+(n_async_tasks*2);
    } else{
        total_reports=5;
    }
    if(verbose){
        if (!(verbose)(1, total_reports))
            return false;
    }
    bool good=true;
    if (n_async_tasks>1) {
        const int work= self_rows / n_async_tasks;
        std::vector<OtherSearchAsyncTask >tasks;
        std::vector<std::future<void>> futures;
        int reports=1;
        FncProgress fp=nullptr;
        if (verbose){
            fp=[&](const int iter, const int n_iter){
                bool continuing=(verbose)(++reports, total_reports);
                return continuing;
            };
        }
        for (int i=0;i<n_async_tasks;i++) {
            const int start=i*work, end = i<n_async_tasks-1?start+work:self_rows;
            tasks.push_back(OtherSearchAsyncTask(
                    nn_graph, nn_graph_columns, otherData, indptr, indices,
                    selfData, self_rows, self_columns, dist, fp, start, end));
        }
        for (OtherSearchAsyncTask &task:tasks) {
            auto f = std::async(task);
            futures.push_back(std::move(f));
        }
        int goods=0;
        for (auto &f:futures){
            f.get();
        }
    }else{
        good=KnnDescent::initialized_nnd_search(
                nn_graph, nn_graph_columns, otherData, indptr, indices,
                selfData, self_rows, self_columns, dist, verbose, -1, self_rows);
    }
    KnnGraph::sort(nn_graph, self_rows, nn_graph_columns);
    return good;
}

bool KnnDescent::initialized_nnd_search(
        double ***const nn_graph, const int nn_graph_columns,
        double **const other_data, const int *const indptr, const int *const indices,
        double **const self_data, const int self_rows, const int self_columns,
        DistancePtr dist, FncProgress verbose, const int start, const int end) {
    int begin, reportModulus, max_loop_reports;
    if (start==-1){
        begin=0;
        max_loop_reports=3;
        reportModulus=end/4;
    } else {
        begin=start;
        max_loop_reports=1;
        reportModulus=(end-start)/2;
    }
    int reports=0;
    for(int i = begin; i < end; i++){
        std::unordered_set<int> tried;
        for(int col = 0; col < nn_graph_columns; col++){
            tried.insert(int(nn_graph[0][i][col]));
        }
        if(verbose){
            if ((i+1)%reportModulus==0) {
                if (reports < max_loop_reports) {
                    reports++;
                    if (!(verbose)(reports+1, max_loop_reports + 2))
                        return false;
                }
            }
        }

        while(true){
            const int vertex = KnnGraph::smallest_flagged(nn_graph, i, nn_graph_columns);
            if(vertex == -1){
                break;
            }
            const int sz=indptr[vertex +1] - indptr[vertex];
            int * candidates= nullptr;
            if (sz>0) {
                candidates=new int[sz];
            } else {
                std::cerr << "Bad input parameters:  "
                             "indptr not ascending from " << vertex << " to " << (vertex+1) << std::endl;
                std::cerr  <<  "Premature termination of initialized_nnd_search... sigh" << std::endl;
                return false;
            }
            int *to=candidates;
            for (const int *start = indices + indptr[vertex], *end = start + sz;
                 start < end; start++, to++) {
                *to = *start;
            }
            for(int j = 0; j < sz;j++){
                if(candidates[j] == vertex || candidates[j] == -1 || tried.find(candidates[j]) != tried.end()){
                    continue;
                }
                const double  d = dist->compute_sortable(other_data[candidates[j]], self_data[i], self_columns);
                numDistances += 1;
                KnnGraph::push(nn_graph, i, d, candidates[j], 1, nn_graph_columns);
                numUncheckedHeapPush += 1;
                tried.insert(candidates[j]);
            }
            delete []candidates;
        }
    }
    if(verbose){
        if (!(verbose)(max_loop_reports + 2, max_loop_reports + 2))
            return false;
    }
    return true;
}
KnnDescent::OtherSearchAsyncTask::OtherSearchAsyncTask(
        double ***const nn_graph, const int nn_graph_columns,
        double **const otherData, const int *const indptr, const int *const indices,
        double **const selfData, const int selfRows, const int selfColumns,
        DistancePtr dist, FncProgress  verbose, const int start, const int end)
        : nn_graph(nn_graph), nn_graph_columns(nn_graph_columns), otherData(otherData), indptr(indptr),
          indices(indices), selfData(selfData), selfRows(selfRows), selfColumns(selfColumns),
          dist(dist), verbose(verbose), start(start), end(end) {
}

void KnnDescent::OtherSearchAsyncTask::operator ()(){
    KnnDescent::initialized_nnd_search(
            nn_graph, nn_graph_columns, otherData, indptr, indices,
            selfData, selfRows, selfColumns, dist, verbose, start, end);
}
KnnDescent::SelfSearchAsyncTask::SelfSearchAsyncTask(
        double ***const nnGraph, const double **const data, const int rows, const int columns,
        const int nNeighbors, long *const rngState, const int maxCandidates, const int nIters,
        const double delta, const double rho, DistancePtr dist,
        double ***old_candidates,  double *** new_candidates,
        const int start, const int end)
        : nn_graph(nnGraph),
          data(data), rows(rows),
          columns(columns),
          n_neighbors(nNeighbors),
          rng_state(rngState),
          max_candidates(
                  maxCandidates),
          n_iters(nIters),
          delta(delta), rho(rho),
          dist(dist),
          old_candidates(old_candidates),
          new_candidates(new_candidates),
          start(start), end(end) {
}

int KnnDescent::SelfSearchAsyncTask::operator()() {
    //nn_descent_internal_low_memory(nn_graph, (const double **)data, rows, columns, n_neighbors, rng_state, max_candidates, n_iters, delta, rho, verbose, start, end);
    int p;
    int q;
    long double d;
    int c = 0;
    int push_in = 1;

    for(int i = start; i < end; i ++) {
        for (int j = 0; j < max_candidates; j++) {
            p = new_candidates[0][i][j];
            if (p < 0) {
                continue;
            }
            for(int k = j; k < max_candidates; k++){
                q = new_candidates[0][i][k];
                if(q < 0){
                    continue;
                }
                const double  d = dist->compute_sortable(data[p], data[q], columns);
                numDistances += 1;
                c += KnnGraph::push(nn_graph, p, d, q, push_in, n_neighbors);
                numHeapPush += 1;
                if(p != q){
                    c += KnnGraph::push(nn_graph, q, d, p, push_in, n_neighbors);
                    numHeapPush += 1;
                }
            }
            for(int k = 0; k < max_candidates; k++){
                q = old_candidates[0][i][k];
                if(q < 0){
                    continue;
                }
                const double  d = dist->compute_sortable(data[p], data[q], columns);
                numDistances += 1;
                c += KnnGraph::push(nn_graph, p, d, q, push_in, n_neighbors);
                numHeapPush += 1;
                if(p != q){
                    c += KnnGraph::push(nn_graph, q, d, p, push_in, n_neighbors);
                    numHeapPush += 1;
                }
            }
        }
    }
    return c;
}
