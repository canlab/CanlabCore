/*
AUTHORS
    Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
    Connor Meehan <connor.gw.meehan@gmail.com>
    Stephen Meehan <swmeehan@stanford.edu>

Provided by suh ( Stanford University's Herzenberg Lab)
License: BSD 3 clause
*/
#include <iostream>
#include <vector>
#include "KnnGraph.h"
#include "KnnDescent.h"
#include "RpTree.h"
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <algorithm>
#include "distances.h"
#include "KnnSearch.h"
//#include "mex.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include<stdio.h>
#include<string.h>
#include <exception>
#include "string.h"

using namespace std::chrono;


int main(int argc, char **argv) {
    std::cout << "Usage:  " << argv[0] << " <n_async_tasks> " << std::endl;
    std::cout << "        optional args 2 & 3: <test_set> <training_set> in ~/Documents/run_umap_examples" << std::endl;
    const bool knn_search_of_self_in_other= argc > 3;
    const int do_async = argc < 2 ? 3 : atoi(argv[1]);
    std::string test_set;

    if (knn_search_of_self_in_other) {
        std::string training_set;
        const bool has_run_umap_examples = argc > 3 && suh::strcmpi(argv[3], "*") != 0;
        if (has_run_umap_examples) {
            test_set = suh::run_umap_examples_file(argv[2]);
            training_set = suh::run_umap_examples_file(argv[3]);
        } else {
            test_set = "../../Data/testSet";
            training_set = "../../Data/trainingSet";
        }
        // The TWO csv files loaded are the 1st and smallest test possible (in Python) for a nn_descent style
        // knnsearch of the nearest neighbors of one data set (self) in another (other)
        // Our C++ wrapper class KnnSearch is referred to as the algorithm for nnd on self and other.
        //
        // In Connor's code for UMAP.m the matlab simple equivalent of all this work can be seen
        // on line #644:
        //[indices, dists] = knnsearch(U.raw_data,X,'K',U.n_neighbors,'Distance',U.metric);

        // The included data for this test is in the nndescent/Data folder:
        //      1. trainingSet.csv (other data matrix without training label @ end column)
        //      2. trainingSet.indptr.csv (self_graph indptr computed on trainingSet self search results)
        //      3. trainingSet.indices.csv (self_graph indices computes on trainingSet self search results)
        //      4. testSet.csv (self data) (other data matrix without training label @ end column)
        //      5. testSetKnnIndices.csv (correct indices answer for nearest neighbors in training set for test set)
        //      6. testSetKnnDists (correct distances answer for nearest neighbors in training set for test set)
        //      5. testSetKnnIndices.csv (correct indices answer for nearest neighbors in training set for test set)
        //      6. testSetKnnDists (correct distances answer for nearest neighbors in training set for test set)
        //
        //Additional test cases of datasets can be generated in MatLab with NnDescent.MakeTestExample(training_set_file, test_set_file)
        //  the 2 parameters specify the location of the csv files without the csv extension.  From these arguments
        //  NnDescent.MakeTestExample
        //      1. expects the two csv files to include the training set label on the last column
        //      2. generates all equivalent files to those described above for the 1st and smallest test
        //
        //To run Python on this exact same data
        //      1. Use the script nndescent/Python/doUmap.py
        //      2. Set current directory for script nndescent/Python (our prefix path below)
        //      3. first CREATE the umap template
        //          doUmap.py trainingSet.csv --n_neighbors=15 --metric=euclidean --min_dist 0.3 --firstRow 0 --verbose --saveTemplate trainingSet.template --output_dimensions 2 --labels trainingSet.labels.csv
        //      4.  APPLY the umap template to run initialise_search and initialized_nn_descent with
        //          doUmap.py testSet.csv --n_neighbors 15 --metric=euclidean  --min_dist 0.3 --firstRow 0 --verbose  --output_dimensions 2 --useTemplate trainingSet.template.umap


        suh::KnnSelfOtherUseCaseCsv testCase1(training_set,test_set);
        suh::KnnSearch h;
        //h->inputDataCheck();
        h.setSelf(testCase1.self.shared_ptr());
        h.setOther(
                testCase1.other.shared_ptr(),
                testCase1.indptr.shared_ptr(),
                testCase1.indices.shared_ptr());
        h.setCorrectAnswer(
                testCase1.correctKnnIndices.shared_ptr(),
                testCase1.correctKnnDists.shared_ptr());
        h.selectNumNeighbors(15);
        h.selectAlgo(suh::KnnSearch::Algos::NndSelfOther);
        int called=0;
        h.setVerbose([&](const int n, const int n_iters){
            std::cout << "\t" << n << " / " << n_iters << std::endl;
            called++;
            return true;
        });

        h.execute(do_async);
        std::cout << called << " \"verbose\" callbacks called "<< std::endl;

        h.reportTotalTime();
        h.printCorrectAnswerMatching();
    } else {
        const bool has_run_umap_examples = argc > 2 && suh::strcmpi(argv[2], "*") != 0;
        if (has_run_umap_examples) {
            test_set = suh::run_umap_examples_file(argv[2]);
        } else {
            test_set =  "../../Data/trainingSet";
        }
        int called=0;
        suh::CsvMatrix<double> testSet(test_set+".csv"),
                correctKnnIndices(test_set+".knnIndices.csv"),
                correctKnnDists(test_set+".knnDists.csv");
        suh::KnnSearch h;
        h.setVerbose([&](const int n, const int n_iters){
            std::cout << "\t" << n << " / " << n_iters << std::endl;
            called++;
            return true;
        });
        h.setSelf(testSet.shared_ptr());
        h.setCorrectAnswer(correctKnnIndices.shared_ptr(), correctKnnDists.shared_ptr());
        h.selectNumNeighbors(15);
        h.selectMetric("Mahalanobis");
        h.selectAlgo(suh::KnnSearch::Algos::NndSelf);
        h.execute(do_async);
        std::cout << called << " \"verbose\" callbacks called "<< std::endl;
        h.reportTotalTime();
        h.printCorrectAnswerMatching();
        if (argc > 3) {
            h.selectCollectStats(true);

            for (int x = 0; x < 2; x++) {
                if (x == 0) {
                    h.selectLowMemory(true);
                } else {
                    h.selectLowMemory(false);
                }
                for (int y = 1; y < 3; y++) {
                    h.selectNumNeighbors(y * 15);
                    for (int j = 1; j < 3; j++) {
                        h.selectMax_candidates(25 * j);
                        h.execute();
                    }
                }
            }
        }
    }
    exit(0);
}
