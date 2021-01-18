/*
AUTHORSHIP
    Connor Meehan <connor.gw.meehan@gmail.com>
    Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
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
#include "mex.h"
#include <iostream>
#include <vector>
#include "KnnGraph.h"
#include "KnnDescent.h"
#include "RpTree.h"
#include "distances.h"
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <algorithm>
#include "distances.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include<stdio.h>
#include<string.h>
#include <exception>

//#define DEBUG 1

using namespace std::chrono;
using namespace suh;

std::string mxGetClassName(const mxClassID classId){
    mxArray *m=mxCreateNumericMatrix(1, 1, classId, mxREAL);
    std::string s(mxGetClassName(m));
    mxDestroyArray(m);
    return s;
}

std::string strNumericType(const mxArray *arg, const bool isVector){
    std::string type;
    if (isVector){
        type="1D ";
        type += mxGetClassName(arg);
        type += " vector";
    } else {
        type="2D ";
        type += mxGetClassName(arg);
        type += " matrix";
    }
    return type;
}

std::vector<int> argNumericSize(const mxArray *arg, const char*name, 
        const bool isVector, const int min_length,  
        const int max_length=2000000000, 
        const mxClassID classId=mxDOUBLE_CLASS){
    int M=-1, N=-1;
    if (!mxIsNumeric(arg) || classId!=mxGetClassID(arg)){
        std::cerr << "parameter " << name << " is a " 
                << strNumericType(arg, isVector)
                << " but needs to be " << mxGetClassName(classId) << ", ";
    } else {
        int M_ = mxGetM(arg), N_ = mxGetN(arg);
        int length=M_*N_;
        if (length<min_length || length > max_length){
            if (min_length != max_length){
                std::cerr <<"length("  << name << "[" << M_ <<"]["
                        << N_ <<"]) must be >=" << min_length
                        << "&&<="<< max_length  << ", ";
            } else if (min_length ==1){
                std::cerr << "parameter "  << name << "[" << M_ <<"]["
                        << N_ <<"]) must be a scalar, ";
            } else {
                std::cerr<<"length("  << name << "[" << M_ <<"]["
                        << N_ <<"]) must be " << min_length << ", ";
            }
        } else {
            if (!isVector || M_<=1 || N_<=1){
                M=M_;
                N=N_;
            }else{
                std::cerr << name << "[" << M_ <<"][" << N_ <<"] " 
                    << strNumericType(arg, isVector)  
                    << " must be a vector, ";
            }
        }
    }
    std::vector<int>sz;
    sz.push_back(M);
    sz.push_back(N);
    return sz;
}


bool isArgFunction(const mxArray *arg, const char *name){
    if (mxFUNCTION_CLASS != mxGetClassID(arg)){
        std::cerr << "parameter " << name << " is a " 
                << mxGetClassName(arg) << " but needs to be " 
                << mxGetClassName(mxFUNCTION_CLASS) << ", ";
        return false;
    }
    return true;
}

bool isArgScalarInRange(const mxArray *arg, const char*name, 
        const double min,  const double max, 
        const mxClassID classId=mxDOUBLE_CLASS){
    std::vector<int> sz=argNumericSize(arg, name, true, 1, 1, classId);
    if (sz[0]==1 && sz[1]==1){
        double v=*mxGetPr(arg);
        if (v<min || v>max){
            std::cerr<< "parameter " << name  << 
                    " must be>=" << min<< " and <=" << max <<", ";
            return false;
        }
        return true;
    }
    return false;
}

int getDistanceBadArgs(const mxArray *prhs[]);
const KnnGraphPtr execute(MatrixPtr  self, const int n_neighbors, 
        const char*distance_metric,  void *distance_args,
        const int n_async_tasks,  const FncProgress fnc_progress);

const KnnGraphPtr execute(MatrixPtr self, MatrixPtr other,
        MatrixIntPtr indptr, MatrixIntPtr indices, const int n_neighbors, 
        const char *distance_metric,  void *distance_args, int n_async_tasks, 
        double transform_queue_size, const FncProgress fnc_progress);

bool progress( mxArray **args, const int iter, const int n_iter){
    int *progress_ptr=(int*)mxGetData(args[1]);
    progress_ptr[0]=iter;
    progress_ptr[1]=n_iter;
    mxArray *argout[1];
    if (mexCallMATLABWithTrap(1, argout, 2, args, "feval") == NULL){
        int *continuing=(int*)mxGetPr(argout[0]);
        return continuing[0]==1;
    }
    mexPrintf("Callback exception occurred... but continuing....\n");
    return true;
}

void showUsageAndExit(){
    const char *core_args="nn_descent self_data, n_neighbors, distance_metric, distance_args, rand_seed, n_async_tasks, ";
    std::cerr << "Self search usage:  " << core_args 
            << "progress_callback" << std::endl;
    std::cerr << "Self/other usage:  " << core_args 
            << "other_data, " << std::endl << "\t\tother_indptr, other_indices, "
            << "transform_queue_size, progress_callback" << std::endl;
    std::cerr << std::endl << "Note:  progress_callback is OPTIONAL"
            << std::endl;
    mexErrMsgTxt("Incorrect input parameter(s)!");
}


#define B_OUT1 plhs[0]
#define B_OUT2 plhs[1]

// nn search within self
#define SELF_IN prhs[0]
#define NN_IN prhs[1]
#define DIST_METRIC_IN prhs[2]
#define DIST_ARGS_IN prhs[3]
#define RAND_IN prhs[4]
#define ASYNC_TASKS_IN prhs[5]
#define FNC_SELF_IN prhs[6]

// nn search of self with other
#define OTHER_IN prhs[6]
#define INDPTR_IN prhs[7]
#define INDICES_IN prhs[8]
#define TRANSFORM_QUEUE_SIZE_IN prhs[9]
#define FNC_OTHER_IN prhs[10]


void checkArgs(const int nrhs, const mxArray *prhs[]){
    /* Check the number of parameters */
    if (nrhs != 6 && nrhs != 7 && nrhs != 10 && nrhs != 11 ) {
        showUsageAndExit();
    } 
    int bad_args=0;
    
#define CHECK_SIZE(Z) if (Z[0]<0|| Z[1] < 0)bad_args++;
    
    std::vector<int>self_size=argNumericSize(SELF_IN, "self_data", false,
            10);
    CHECK_SIZE(self_size);
    if (!isArgScalarInRange(NN_IN, "n_neighbors", 1.0, 100.0))
        bad_args++;
    bad_args+=getDistanceBadArgs(prhs);
    if (!isArgScalarInRange(RAND_IN, "rand_seed", -1000, 1000.0))bad_args++;
    if (!isArgScalarInRange(ASYNC_TASKS_IN, "n_async_tasks", 1.0, 100.0))bad_args++;
    if (nrhs<8){
        if (nrhs==7 && !isArgFunction(FNC_SELF_IN, "progress_callback" )) 
            bad_args++;
    }else{
        std::vector<int>other_size=argNumericSize(OTHER_IN, 
                "other_data", false, 10);
        CHECK_SIZE(other_size);
        if (self_size[1] != other_size[1]){
            bad_args++;
            #define OUT_2D(Z) "[" << Z[0]<<  "][" << Z[1] << "]" << (Z[0] * Z[1])

            std::cerr << "columns mismatch self_data=" << OUT_2D(self_size)
                    << " but other_data" << OUT_2D(other_size) << ", ";
        }
        const int sz=other_size[0];
        std::vector<int>indptr_size=argNumericSize(INDPTR_IN, 
                "other_indptr", true, sz+1, sz+1, mxINT32_CLASS);
        CHECK_SIZE(indptr_size);
        int n_neighbors = (int)*mxGetPr(NN_IN);
        std::vector<int>indices_size=argNumericSize(INDICES_IN,
                "other_indices", true, sz+1, n_neighbors*n_neighbors*(sz+1), 
                mxINT32_CLASS);
        CHECK_SIZE(indices_size);
        if (!isArgScalarInRange(TRANSFORM_QUEUE_SIZE_IN, 
                "transform_queue_size", 1.0, 4.0))
            bad_args++;
        if (nrhs==11 && !isArgFunction(FNC_OTHER_IN, "progress_callback")) bad_args++;
    }
    if (bad_args>0){
        std::cerr<<std::endl<<bad_args<<" bad parameter(s) "<<std::endl;
        showUsageAndExit();
    }
}

bool isEmpty(const mxArray *arg){
    const int M=mxGetM(arg), N=mxGetN(arg);
    return M< 1 || N< 1;
}

MatrixPtr getMatrix(const mxArray *arg){
    double *ptr=mxGetPr(arg);
    const int M=mxGetM(arg), N=mxGetN(arg);
    Matrix<double>matrix(ptr, M, N, suh::Transfer::kCopyColumnWise);
    return matrix.shared_ptr();
}

bool isScalar(const mxArray *arg){
    const int M= mxGetM(arg), N=mxGetN(arg);
    return M==1 && N==1;
}

#define IS(r) suh::strcmpi(distance_metric,r)==0

int getDistanceBadArgs(const mxArray *prhs[]){
    int bad_args=0;
    char *distance_metric=nullptr;
    if (!mxIsChar(DIST_METRIC_IN)){
        bad_args++;
        std::cerr<< "distance_metric must be char array";
    }else{
        distance_metric=mxArrayToString(DIST_METRIC_IN);
        if (!DistanceMetric::IsSupported(distance_metric)){
            bad_args++;
            std::cerr << "distance_metric \"" << distance_metric
                    <<"\" not supported, ";
        }
    }
    if (distance_metric != nullptr){        
        std::vector<int>Z=argNumericSize(DIST_ARGS_IN,
                "distance_args", false, 0);
        if (Z[0]<0|| Z[1] < 0){
            bad_args++;
        }else if (Z[0]>0 && Z[1] >0) {            
            if (IS("minkowski")){  
                if (Z[0] != 1 || Z[1] != 1){
                    std::cerr << "Minkowski distance_args size is "
                            << Z[0] << " X " << Z[1]
                            << " but must be a scalar, ";
                    bad_args++;
                } else {
                    const double * distance_args_ptr = mxGetPr(DIST_ARGS_IN);
                    if (*distance_args_ptr<=1.0){
                        std::cerr << "Minkowski distance_args " 
                                << *distance_args_ptr << 
                                " must >= 1.0, ";
                        bad_args++;
                    }
                }
            }else if (IS("mahalanobis")){
               if (Z[0] != Z[1]){
                    std::cerr << "Mahalanobis distance_args size is "
                            << Z[0] << " X " << Z[1]
                            << " but must be a SQUARE covaraiance "
                            << "matrix for self_data, ";
                    bad_args++;
                } else {
                    const int self_M = mxGetM(SELF_IN), self_N = mxGetN(SELF_IN);
                    if (self_N != Z[1]){
                        std::cerr <<"Mahalanobis distance_args size is "
                                << Z[0] << " X " << Z[1] << " but must be "
                                << self_N  << " X " <<self_N 
                                << " since self_data is " 
                                << self_M  << " X "<<self_N;
                        bad_args++;
                    }
                }
            }else if (IS("seuclidean")){
               if (Z[0] != 1 ){
                    std::cerr << "SEuclidean distance_args size is "
                            << Z[0] << " X " << Z[1]
                            << " but must be a single row"
                            << "of deviations from self_data, ";
                    bad_args++;
                } else {
                    const int self_M = mxGetM(SELF_IN), self_N = mxGetN(SELF_IN);
                    if (self_N != Z[1]){
                        std::cerr <<"SEuclidean distance_args size is "
                                << Z[0] << " X " << Z[1] 
                                << " but must be 1 X " << self_N 
                                << " since self_data is " 
                                << self_M  << " X "<<self_N;
                        bad_args++;
                    }
                }
            } else {
                bad_args++;
                std::cerr << "distance_args should be empty for  \"" 
                        << distance_metric <<"\", ";
            }
        }
        //mxFree(distance_metric);
    }
    return bad_args;
}
    
void mexFunction(
        int nlhs, mxArray *plhs[], /* Output variables */
        int nrhs, const mxArray *prhs[]) /* Input variables */
{
    checkArgs(nrhs, prhs);
    const int self_M = mxGetM(SELF_IN), self_N = mxGetN(SELF_IN);
    const double * self_ptr = mxGetPr(SELF_IN);
    int n_neighbors = (int)*mxGetPr(NN_IN);
    char *distance_metric=mxArrayToString(DIST_METRIC_IN);
    void *distance_args_ptr = nullptr;
    MatrixPtr matrix_dist_args;
    if (!isEmpty(DIST_ARGS_IN)){
        if (IS("minkowski")){
            distance_args_ptr=mxGetPr(DIST_ARGS_IN);
        } else if (IS("mahalanobis")){
            matrix_dist_args=getMatrix(DIST_ARGS_IN);
            distance_args_ptr=matrix_dist_args->matrix_;
        } else if (IS("seuclidean")){
            matrix_dist_args=getMatrix(DIST_ARGS_IN);
            distance_args_ptr=matrix_dist_args->matrix_;
        }
    }
    
    int rand_seed=(int) mxGetScalar(RAND_IN);
    if (rand_seed == 0) {
        srand(0);
    }
    int n_async_tasks= (int)*mxGetPr(ASYNC_TASKS_IN);
#ifdef DEBUG
    std::cout << "n_neighbors=" << n_neighbors << ", distance_metric=\""
            << distance_metric << "\", rand_seed=" << rand_seed
            << ", n_async_tasks="
            << n_async_tasks << std::endl;
#endif

    FncProgress fnc_progress=nullptr;
    mxArray *args[2];
    Matrix<double>self(self_ptr, self_M, self_N, suh::Transfer::kCopyColumnWise);
    KnnGraphPtr nn_graph;
    if (nrhs>7){
        if (nrhs==11){
            args[0]=(mxArray *) FNC_OTHER_IN;
            args[1]=mxCreateNumericMatrix(1, 2, mxINT32_CLASS, mxREAL);
            fnc_progress=[&](const int iter, const int n_iter){
                bool continuing=progress(args, iter, n_iter);
                return continuing;
            };
        }

        const double * other_ptr = mxGetPr(OTHER_IN);
        const int other_M = mxGetM(OTHER_IN), other_N = mxGetN(OTHER_IN);
        Matrix<double>other(other_ptr, other_M, other_N,
                            suh::Transfer::kCopyColumnWise);

        const int indptr_M=mxGetM(INDPTR_IN), indptr_N=mxGetN(INDPTR_IN);
        const int *indptr_ptr=(int *)mxGetData(INDPTR_IN);
        Matrix<int> indptr(indptr_ptr, 1, indptr_M, suh::Transfer::kCopy);


        const int indices_M=mxGetM(INDICES_IN), indices_N=mxGetN(INDICES_IN);
        const int *indices_ptr=(int *)mxGetData(INDICES_IN);
        Matrix<int> indices(indices_ptr, 1, indices_M, suh::Transfer::kCopy);

        const double *transform_queue_size=mxGetPr(TRANSFORM_QUEUE_SIZE_IN);
        nn_graph=execute(
                self.shared_ptr(), other.shared_ptr(),
                indptr.shared_ptr(), indices.shared_ptr(),
                n_neighbors, distance_metric, distance_args_ptr,
                n_async_tasks, *transform_queue_size, fnc_progress);
    } else{
        if (nrhs==7){
            args[0]=(mxArray *) FNC_SELF_IN;
            args[1]=mxCreateNumericMatrix(1, 2, mxINT32_CLASS, mxREAL);
            fnc_progress=[&](const int iter, const int n_iter){
                bool continuing=progress(args, iter, n_iter);
                return continuing;
            };
        }
        nn_graph=execute(self.shared_ptr(), n_neighbors, 
                distance_metric, distance_args_ptr, 
                n_async_tasks, fnc_progress);
    }
    if (nn_graph){
        B_OUT1 = mxCreateDoubleMatrix(self_M, n_neighbors,mxREAL);
        B_OUT2 = mxCreateDoubleMatrix(self_M, n_neighbors,mxREAL);
        nn_graph->copy_indices_MatLab(mxGetPr(B_OUT1));
        nn_graph->copy_distances_column_wise(mxGetPr(B_OUT2));
    } else{
        mexPrintf("No knn graph ... progress callback halted !!");
        B_OUT1 = mxCreateDoubleMatrix(0, 0,mxREAL);
        B_OUT2 = mxCreateDoubleMatrix(0, 0,mxREAL);
    }
    mxFree(distance_metric);
}


const KnnGraphPtr execute(
        MatrixPtr self, MatrixPtr other, MatrixIntPtr indptr, 
        MatrixIntPtr indices, const int n_neighbors, 
        const char*distance_metric, void *distance_args,
        int n_async_tasks, double transform_queue_size, 
        FncProgress fnc_progress){
    long rng_state[3];
    rng_state[0] = rand();
    rng_state[1] = rand();
    rng_state[2] = rand();
    if (fnc_progress==nullptr){
        #ifdef DEBUG
        mexPrintf("self in other n_async_tasks=%d\n", n_async_tasks);
        #endif
    }else
        n_async_tasks=1;
    
    if (transform_queue_size > 4.0)
        transform_queue_size=4.0;
    else if (transform_queue_size<1.0)
        transform_queue_size=1.0;
    #ifdef DEBUG
    mexPrintf("self in other transform_queue_size=%f\n", transform_queue_size);
    #endif
    return KnnDescent::search(
            other, self, indptr, indices, n_neighbors, 
            (float)transform_queue_size,  rng_state, distance_metric, 
            n_async_tasks, fnc_progress, distance_args);
}

const KnnGraphPtr execute(MatrixPtr self, const int n_neighbors,
        const char*distance_metric, void *distance_args, 
        const int n_async_tasks, const FncProgress fnc_progress){
    long rng_state[3];
    rng_state[0] = rand();
    rng_state[1] = rand();
    rng_state[2] = rand();
#ifdef DEBUG
    if (fnc_progress==nullptr){
        mexPrintf("self n_async_tasks=%d\n", n_async_tasks);
    }
#endif
    return KnnDescent::search(self,n_neighbors,rng_state,
            distance_metric, n_async_tasks, fnc_progress, 
            distance_args);
}
