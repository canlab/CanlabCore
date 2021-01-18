//
// Created by Stephen Meehan on 10/29/20.
//
/*
AUTHORS
    Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
    Connor Meehan <connor.gw.meehan@gmail.com>
    Stephen Meehan <swmeehan@stanford.edu>

Provided by suh ( Stanford University's Herzenberg Lab)
License: BSD 3 clause
*/

#include "KnnSearch.h"
#include "KnnDescent.h"
#include "KnnGraph.h"
#include "RpTree.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <chrono>

using namespace std::chrono;

suh::KnnSearch::KnnSearch() {
    /*
    class comma_numpunct : public std::numpunct<char>
    {
    public:
        comma_numpunct(char thousands_sep, const char* grouping)
                :m_thousands_sep(thousands_sep),
                 m_grouping(grouping){}
    protected:
        char do_thousands_sep() const{return m_thousands_sep;}
        std::string do_grouping() const {return m_grouping;}
    private:
        char m_thousands_sep;
        std::string m_grouping;

    };

    std::locale comma_locale(std::locale(), new comma_numpunct(',', "\03"));
    std::cout.imbue(comma_locale);*/

    //initializing setup

    /*rng_state=new long[3];
    rng_state[0] = random();
    rng_state[1] = random();
    rng_state[2] = random();*/

    //for exact python testing
    rng_state = new long[3]{394951199, 1210874715, 266162594};

    angular = false;
    n_trees = 0;
    n_neighbors = 15;
    leaf = true;

    max_candidates = 50;
    n_iters = 10;
    delta = 0.001;
    rho = 0.5;
    low_memory = true;
    verbose = nullptr;
    fullStats = false;
    headers = "Algorithm Name,Rows,Cols, Neighbors, Candidates, Iters, Leaf Execution Time, Algo Execution Time, Total Execution Time, Distances, Heap Push, Unchecked Heap Push";
}

void suh::KnnSearch::selectAlgo(int input) {
    if (input == Algos::NndSelf) {
        algoSelected = Algos::NndSelf;
        nameOfAlgo = "nndescent";
    } else if (input == Algos::NndSelfOther) {
        algoSelected = Algos::NndSelfOther;
        nameOfAlgo = "initialized";
    } else {
        std::cout << "Invalid selection" << std::endl;
        //exit(0);
    }
}

void suh::KnnSearch::selectAlgo(std::string input) {
    std::transform(input.begin(), input.end(), input.begin(), ::tolower);
    if (input == "nndescent") {
        algoSelected = Algos::NndSelf;
        nameOfAlgo = input;
    } else if (input == "initialized") {
        nameOfAlgo = input;
        algoSelected = Algos::NndSelfOther;
    } else {
        std::cout << "Invalid selection" << std::endl;
        //exit(0);
    }
}

void suh::KnnSearch::selectMetric(std::string input, const void *dist_arguments) {
    std::transform(input.begin(), input.end(), input.begin(), ::tolower);
    dist_metric = input;
    dist_args = (void *)dist_arguments;
}

void suh::KnnSearch::selectNumNeighbors(int n_neighbors_in) {
    if (n_neighbors_in > 30) {
        std::cout << "This algorithm may not be optimal for the number of neighbors requested" << std::endl;
    }
    n_neighbors = n_neighbors_in;
}

void suh::KnnSearch::selectLeaf(bool leafIn) {
    leaf = leafIn;
}

void suh::KnnSearch::reportTotalTime(){
    if (knn_graph) {
        std::cout << suh::string_format(
                "Pure search %4.3f seconds (total search cost %4.3f  minus RpForest cost %4.3f) ",
                suh::secs(durationAlgo.count() - knn_graph->rp_forest_->get_duration()),
                suh::secs(durationAlgo.count()),
                suh::secs(knn_graph->rp_forest_->get_duration()))
                  << std::endl;
    } else{
        std::cout << "Progress cancelled ... sigh" <<std::endl;
    }
}

void suh::KnnSearch::stats() {
    std::cout << "Algo: " << nameOfAlgo << std::endl;
    std::cout << "Data size: " << selfData->rows_ << "x" << selfData->columns_ << std::endl;
    std::cout << "Neighbors: " << n_neighbors << std::endl;
    if (fullStats) {
        std::cout << "Leaf: " << leaf << std::endl;
        std::cout << "n_iters: " << n_iters << std::endl;
        std::cout << "n_trees: " << n_trees << std::endl;
        std::cout << "max_candidates: " << max_candidates << std::endl;
        std::cout << "delta: " << delta << std::endl;
        std::cout << "rho: " << rho << std::endl;
        std::cout << "low_memory: " << low_memory << std::endl;
        std::cout << "angular: " << angular << std::endl;
    }
    std::cout << "Time taken to build random projection forest: " << durationLeaf << " microseconds" << std::endl;
    std::cout << "Time taken to execute selected algo: " << durationAlgo.count() << " microseconds" << std::endl;
    std::cout << "Total time: " << (durationAlgo.count() + durationLeaf) << " microseconds" << std::endl;
}


void suh::KnnSearch::execute(const int n_async_tasks) {

    bool goodToGo = false;
    if (selfData->rows_ < 2)
        std::cerr << "NOTHING to do self data rows are <2!" << std::endl;
    else if (selfData->columns_ == 0)
        std::cerr << "NOTHING to do self data columns are 0!" << std::endl;
    else
        goodToGo = true;

    cleanStats();


    auto startAlgo = high_resolution_clock::now();
    if (goodToGo) {
        if (algoSelected == Algos::NndSelf) {
            std::cout << "Calling Jonathan's translation from Python to C++ of nn_descent" << std::endl;
            //low_memory=false;
            knn_graph = KnnDescent::search(
                    selfData, n_neighbors, rng_state,dist_metric,  n_async_tasks, verbose,
                    dist_args, max_candidates,n_iters, delta, rho, leaf, low_memory,  n_trees, angular);
        } else if (algoSelected == Algos::NndSelfOther) {
            std::cout
                    << "Calling Jonathan's translation from Python to C++ of initialise_search and initialized_nn_descent"
                    << std::endl;
            this->knn_graph= KnnDescent::search(
                    otherData, selfData, otherIndptr, otherIndices, n_neighbors, transform_queue_size,
                    rng_state, dist_metric,  n_async_tasks,  verbose, dist_args, angular);
        } else if (algoSelected == -1) {
            std::cout << "No algorithm selected" << std::endl;
            // do not want to exit(0) the MatLab application;
        } else {
            std::cout << "Invalid algorithm selected" << std::endl;
        }
    }
    auto stopAlgo = high_resolution_clock::now();
    durationAlgo = duration_cast<microseconds>(stopAlgo - startAlgo);
    if (collectStats) {
        dataCollectionFile << nameOfAlgo << "," << selfData->rows_ << "," << selfData->columns_ << "," << n_neighbors
                           << "," << max_candidates << "," << n_iters << "," << durationLeaf << ","
                           << durationAlgo.count() << "," << (durationAlgo.count() + durationLeaf) << ", "
                           << KnnDescent::numDistances << ", " << KnnDescent::numHeapPush << ", "
                           << KnnDescent::numUncheckedHeapPush << std::endl;
    }
}

void suh::KnnSearch::printCorrectAnswerMatching() {
    if (correctKnnIndices == nullptr || !knn_graph) return;
    suh::MatrixPtr inds = knn_graph->indices();
    if (correctKnnIndices->rows_ != inds->rows_ ||
        correctKnnIndices->columns_ != inds->columns_) {
        if (correctKnnIndices->columns_ > 0)
            std::cerr << "Matrix rows/columns MUST be same size" << std::endl;
        return;
    }
    suh::print_inequalities(*inds, *correctKnnIndices, true);
    suh::MatrixPtr dists = knn_graph->distances();
    suh::print_inequalities(*dists, *correctKnnDists, true);
}

void suh::KnnSearch::selectFullStats(const bool in) {
    fullStats = in;
}

suh::KnnSearch::~KnnSearch() {
//   KnnGraph::Deallocate((double ***)nn_graph);

}

void suh::KnnSearch::printIndices(const int row) const {
    if (!knn_graph) return;
    std::cout << " inidices@#" << row << ": ";
    double ***nn_graph=*this->knn_graph;
    for (int neighbor = 0; neighbor < n_neighbors; neighbor++) {
        const int idx = (int) nn_graph[0][row][neighbor];
        std::cout << idx << " ";
    }
    std::cout << std::endl;
}

void suh::KnnSearch::printDistances(const int row) const {
    if (!knn_graph) return;
    std::cout << "distances@#" << row << ": ";
    double ***nn_graph=*this->knn_graph;
    for (int neighbor = 0; neighbor < n_neighbors; neighbor++) {
        std::cout << nn_graph[1][row][neighbor] << " ";
    }
    std::cout << std::endl;
}


void suh::KnnSearch::inputDataCheck() {
    for (int x = 0; x < selfData->columns_; x++) {
        std::cout << selfData->matrix_[0][x] << std::endl;
    }
}


void suh::KnnSearch::selectCollectStats(bool collectStatsIn) {
    collectStats = collectStatsIn;
    if (collectStats) {
        dataCollectionFile.open("./DataCollection/dataCollection.csv");
        dataCollectionFile << headers << std::endl;
    }
}

void suh::KnnSearch::selectCollectStats(bool collectStatsIn, std::string fileNameIn) {
    collectStats = collectStatsIn;
    if (collectStats) {
        dataCollectionFile.open("./DataCollection/" + fileNameIn);
        dataCollectionFile << headers << std::endl;
    }
}

void suh::KnnSearch::cleanStats() {
    KnnDescent::numDistances = 0;
    KnnDescent::numHeapPush = 0;
    KnnDescent::numUncheckedHeapPush = 0;
}

void suh::KnnSearch::selectN_iters(int itersIn) {
    n_iters = itersIn;
}

void suh::KnnSearch::selectMax_candidates(int candidatesIn) {
    max_candidates = candidatesIn;
}

void suh::KnnSearch::selectLowMemory(bool lowMemoryIn) {
    low_memory = lowMemoryIn;
}



