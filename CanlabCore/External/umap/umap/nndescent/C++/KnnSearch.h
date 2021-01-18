//
// Created by Stephen Meehan on 10/29/20.
//
/*
AUTHORS
   Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
   Stephen Meehan <swmeehan@stanford.edu>
   Connor Meehan <connor.gw.meehan@gmail.com>

Provided by suh ( Stanford University's Herzenberg Lab)
License: BSD 3 clause
*/

#ifndef KNNSEARCH_H
#define KNNSEARCH_H
#include "suh.h"

#include <string>
#include <vector>
#include <chrono>
#include "KnnGraph.h"

using namespace std::chrono;
namespace suh {
    class KnnSearch {
    private:
        std::string dist_metric="euclidean";

        suh::MatrixPtr  selfData, otherData, correctKnnIndices, correctKnnDists;
        suh::MatrixIntPtr otherIndptr, otherIndices;
        float transform_queue_size=1;
        KnnGraphPtr knn_graph;
        std::vector<std::vector<int> > leaf_array;

        //timing
        long durationLeaf;
        microseconds durationAlgo;

        //Collect Data
        std::ofstream dataCollectionFile;
        std::string dataCollectionFileName;
        std::string headers;
        void cleanStats();
        bool collectStats;

        //selected info
        void *dist_args= nullptr;
        std::string nameOfAlgo;
        int algoSelected;
        int max_candidates;
        int n_iters;
        double delta;
        double rho;
        bool low_memory;
        suh::FncProgress verbose;
        long *rng_state;
        bool angular=false;
        int n_trees=-1;
        int n_neighbors;
        bool leaf;
        bool fullStats;

    public:
        inline void setVerbose(FncProgress verbose){
            this->verbose=verbose;
        }
        enum Algos {
            // currently only KnnSearch algorithms are nearest neighbor descent.
            // Other ones will be added in the future to cater to different
            // data patterns with the best balance of accuracy and speed
            NndSelf = 0,
            NndSelfOther,
        };

        KnnSearch();

        void inputDataCheck();

        void selectAlgo(std::string);

        void selectAlgo(int);

        void selectNumNeighbors(int);

        void selectFullStats(bool);

        void selectLeaf(bool);

        void stats();

        void execute(const int n_async_tasks=3);

        void selectCollectStats(bool);

        void selectN_iters(int);

        void selectMax_candidates(int);

        void selectCollectStats(bool,std::string);

        void selectLowMemory(bool);

        void printIndices(const int row) const;
        void printDistances(const int row) const;
        inline void setSelf(std::string csvFile){
            suh::CsvMatrix<double> testSet(csvFile + ".csv");
            setSelf(testSet.shared_ptr());
        }
        inline void setSelf(suh::MatrixPtr selfData){
            this->selfData=selfData;
        }
        inline void setOther(suh::MatrixPtr otherData, suh::MatrixIntPtr otherIndptr,
                             suh::MatrixIntPtr otherIndices){
            this->otherData=otherData;
            this->otherIndptr=otherIndptr;
            this->otherIndices=otherIndices;
        }
        inline void setCorrectAnswer(
                suh::MatrixPtr correctIndices, suh::MatrixPtr correctDists){
            this->correctKnnIndices=correctIndices;
            this->correctKnnDists=correctDists;
        }
        void printCorrectAnswerMatching();
        ~KnnSearch();

        void reportTotalTime();

        void selectMetric(std::string input, const void *dist_args= nullptr);
    };
}

#endif //SUH_KNNSEARCH_H
