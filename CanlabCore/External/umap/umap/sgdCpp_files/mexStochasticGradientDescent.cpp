//
// Created by Connor on 2020-04-28.
//
/*
 *
   AUTHORSHIP
   Primary Developers:  Connor Meehan <connor.gw.meehan@gmail.com>
   			 Stephen Meehan <swmeehan@stanford.edu>
   Math Lead:  		 Connor Meehan <connor.gw.meehan@gmail.com>
   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
   Provided by the Herzenberg Lab at Stanford University
   License: BSD 3 clause
*/


#include "mex.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <exception>

int callProgressFcn( mxArray **args, const int epoch){
    int*epochsPtr=(int*)mxGetData(args[2]);    
    epochsPtr[0]=epoch;
    mxArray *argout[1];
    mexCallMATLAB(1, argout, 3, args, "feval");  /* Call the plotfcn function handle */
    int *continuing=(int*)mxGetPr(argout[0]);
    return continuing[0];
}

void doMoveOther(
        double *head_embedding,  double *tail_embedding,
        const unsigned *head, const unsigned *tail, const int n_epochs, const int n_vertices,
        const double *epochs_per_sample, const double a, const double b,
        const double gamma, const double initial_alpha,
        const double negative_sample_rate, const int n_1_simplices,
        const int n_components, mxArray **prArgs, const int epochReports) {
    double alpha;
    double BG2S;
    double ABNEG2;
    double BNEG1;
    double *epochs_per_negative_sample;
    double *epoch_of_next_negative_sample;
    double *epoch_of_next_sample;
    double nTh;
    int n_epoch;
    alpha = initial_alpha;
    BG2S = 2 * gamma * b;
    ABNEG2 = -2.0 * a * b;
    BNEG1 = b - 1;
    epoch_of_next_negative_sample = static_cast<double *>(mxMalloc(n_1_simplices * sizeof(double)));
    epochs_per_negative_sample = static_cast<double *>(mxMalloc(n_1_simplices * sizeof(double)));
    epoch_of_next_sample = static_cast<double *>(mxMalloc(n_1_simplices * sizeof(double)));

    for (int i = 0; i < n_1_simplices; i++) {
        epochs_per_negative_sample[i] = epochs_per_sample[i] / negative_sample_rate;
    }
    for (int i = 0; i < n_1_simplices; i++) {
        epoch_of_next_negative_sample[i] = epochs_per_negative_sample[i];
        epoch_of_next_sample[i] = epochs_per_sample[i];
    }
    nTh = (double) n_epochs / double(epochReports);
    n_epoch = 1;
    double *current = static_cast<double *>(mxMalloc(n_components * sizeof(double)));
    double *other = static_cast<double *>(mxMalloc(n_components * sizeof(double)));
    double *grad = static_cast<double *>(mxMalloc(n_components * sizeof(double)));
    double *sub = static_cast<double *>(mxMalloc(n_components * sizeof(double)));
    int iRandi = 0;
    int n_neg_samples = 0;
    double grad_coef = 0;
    double dist_squared = 0;
    double val = 0;
    double alpha4 = alpha * 4, alphaNeg4 = alpha * -4;
    for (size_t m = 0; m < n_components; m++) {
        current[m] = 0;
        other[m] = 0;
        grad[m] = 0;
        sub[m] = 0;
    }
    for (int n = n_epoch;n <=n_epochs;n++) {
        for (int i = 0;i < n_1_simplices;i++) {
            if (epoch_of_next_sample[i] > n) {
                continue;
            }
            const int j = head[i] - 1;//const
            int k = tail[i] - 1;
            for (size_t m = 0; m < n_components; m++) {
                current[m] = head_embedding[j + n_vertices*m];
                other[m] = tail_embedding[k + n_vertices*m];
                sub[m] = current[m] - other[m];
            }
            dist_squared = 0;
            for (size_t m = 0; m < n_components; m++) {
                dist_squared += sub[m] * sub[m];
            }
            if (dist_squared > 0) {
                grad_coef = (ABNEG2 * pow(dist_squared, BNEG1)) / (a * pow(dist_squared, b) + 1);
                for (size_t m = 0; m < n_components; m++) {
                    val = grad_coef * sub[m];
                    if (val >= 4) {
                        grad[m] = alpha4;
                    } else if (val <= -4) {
                        grad[m] = alphaNeg4;
                    } else {
                        grad[m] = val * alpha;
                    }
                    current[m] = current[m] + grad[m];
                }
                for (size_t m = 0; m < n_components; m++) {
                    other[m] = other[m] - grad[m];
                    tail_embedding[k + n_vertices*m] = other[m];
                }
            } else {
                for (size_t m = 0; m < n_components; m++) {
                    grad[m] = 0;
                }
            }
            epoch_of_next_sample[i] += epochs_per_sample[i];
            n_neg_samples = static_cast<int>(floor(((static_cast<double>(n))
            - epoch_of_next_negative_sample[i]) /
                    epochs_per_negative_sample[i]));
            
            for (int p = 0; p < n_neg_samples; p++) {
                k = rand() % n_vertices;
                if (j == k) {
                    continue;
                }
                dist_squared = 0;
                for (size_t m = 0; m < n_components; m++) {
                    other[m] = tail_embedding[k + n_vertices*m];
                    sub[m] = current[m] - other[m];
                    dist_squared += sub[m] * sub[m];
                }
                if (dist_squared > 0) {
                    grad_coef = ((BG2S / (0.001 + dist_squared))) / (a * pow(dist_squared, b) + 1);
                    for (size_t m = 0; m < n_components; m++) {
                        val = grad_coef * sub[m];
                        if (val >= 4) {
                            grad[m] = alpha4;
                        } else if (val <= -4) {
                            grad[m] = alphaNeg4;
                        } else {
                            grad[m] = val * alpha;
                        }
                    }
                } else {
                    for (size_t m = 0; m < n_components; m++) {
                        grad[m] = 4;
                    }
                }
                for (size_t m = 0; m < n_components; m++) {
                    current[m] = current[m] + (grad[m]);
                }
            }
            for (size_t m = 0; m < n_components; m++) {
                head_embedding[j + n_vertices*m] = current[m];
            }
            epoch_of_next_negative_sample[i] += n_neg_samples * epochs_per_negative_sample[i];
        }
        alpha = initial_alpha * (1 - static_cast<double>(static_cast<double>(n) / static_cast<double>(n_epochs)));
        alpha4 = alpha * 4;
        alphaNeg4 = alpha * -4;
        if (prArgs!=NULL){
            double nBynTh = floor(fmod((double) n, nTh));
            if (nBynTh == 0) {
                n_epoch = n + 1;
                if (n_epoch < n_epochs) {
                    int continuing=callProgressFcn(prArgs, n_epoch);
                    if (continuing==0){
                        break;
                    }
                }
            }
        }
    }
    mxFree(current);
    mxFree(other);
    mxFree(grad);
    mxFree(sub);
    mxFree(epoch_of_next_negative_sample);
    mxFree(epochs_per_negative_sample);
    mxFree(epoch_of_next_sample);
}

void doNotMoveOther(
        double *head_embedding,  double *tail_embedding,
        const unsigned *head, const unsigned *tail, int n_epochs, const int n_vertices,
        const double *epochs_per_sample, const double a, const double b,
        const double gamma, const double initial_alpha,
        const double negative_sample_rate, const int n_1_simplices,
        const int n_components, const int size_head_embedding,
        mxArray **prArgs, const int epochReports) {
    double alpha;
    double BG2S;
    double ABNEG2;
    double BNEG1;
    double *epochs_per_negative_sample;
    double *epoch_of_next_negative_sample;
    double *epoch_of_next_sample;
    double nTh;
    int n_epoch;
    alpha = initial_alpha;
    BG2S = 2 * gamma * b;
    ABNEG2 = -2.0 * a * b;
    BNEG1 = b - 1;
    epoch_of_next_negative_sample = static_cast<double *>(mxMalloc(n_1_simplices * sizeof(double)));
    epochs_per_negative_sample = static_cast<double *>(mxMalloc(n_1_simplices * sizeof(double)));
    epoch_of_next_sample = static_cast<double *>(mxMalloc(n_1_simplices * sizeof(double)));

    for (int i = 0; i < n_1_simplices; i++) {
        epochs_per_negative_sample[i] = epochs_per_sample[i] / negative_sample_rate;
    }
    for (int i = 0; i < n_1_simplices; i++) {
        epoch_of_next_negative_sample[i] = epochs_per_negative_sample[i];
        epoch_of_next_sample[i] = epochs_per_sample[i];
    }
    nTh = (double) n_epochs / double(epochReports);
    n_epoch = 1;
    double *current = static_cast<double *>(mxMalloc(n_components * sizeof(double)));
    double *other = static_cast<double *>(mxMalloc(n_components * sizeof(double)));
    double *grad = static_cast<double *>(mxMalloc(n_components * sizeof(double)));
    double *sub = static_cast<double *>(mxMalloc(n_components * sizeof(double)));
    int iRandi = 0;
    int n_neg_samples = 0;
    double grad_coef = 0;
    double dist_squared = 0;
    double val = 0;
    double alpha4 = alpha * 4, alphaNeg4 = alpha * -4;
    for (size_t m = 0; m < n_components; m++) {
        current[m] = 0;
        other[m] = 0;
        grad[m] = 0;
        sub[m] = 0;
    }
    for (int n = n_epoch;n <=n_epochs;n++) {
        for (int i = 0;i < n_1_simplices;i++) {
            if (epoch_of_next_sample[i] > n) {
                continue;
            }
            const int j = head[i] - 1;//const
            int k = tail[i] - 1;
            for (size_t m = 0; m < n_components; m++) {
                current[m] = head_embedding[j + size_head_embedding *m];
                other[m] = tail_embedding[k + n_vertices*m];
                sub[m] = current[m] - other[m];
            }
            dist_squared = 0;
            for (size_t m = 0; m < n_components; m++) {
                dist_squared += sub[m] * sub[m];
            }
            if (dist_squared > 0) {
                grad_coef = (ABNEG2 * pow(dist_squared, BNEG1)) / (a * pow(dist_squared, b) + 1);
                for (size_t m = 0; m < n_components; m++) {
                    val = grad_coef * sub[m];
                    if (val >= 4) {
                        grad[m] = alpha4;
                    } else if (val <= -4) {
                        grad[m] = alphaNeg4;
                    } else {
                        grad[m] = val * alpha;
                    }
                    current[m] = current[m] + grad[m];
                }
            } else {
                for (size_t m = 0; m < n_components; m++) {
                    grad[m] = 0;
                }
            }
            epoch_of_next_sample[i] += epochs_per_sample[i];
            n_neg_samples = static_cast<int>(floor(((static_cast<double>(n))
            - epoch_of_next_negative_sample[i]) /
                    epochs_per_negative_sample[i]));
            
            for (int p = 0; p < n_neg_samples; p++) {
                k = rand() % n_vertices;
                dist_squared = 0;
                for (size_t m = 0; m < n_components; m++) {
                    other[m] = tail_embedding[k + n_vertices*m];
                    sub[m] = current[m] - other[m];
                    dist_squared += sub[m] * sub[m];
                }
                if (dist_squared > 0) {
                    grad_coef = ((BG2S / (0.001 + dist_squared))) / (a * pow(dist_squared, b) + 1);
                    for (size_t m = 0; m < n_components; m++) {
                        val = grad_coef * sub[m];
                        if (val >= 4) {
                            grad[m] = alpha4;
                        } else if (val <= -4) {
                            grad[m] = alphaNeg4;
                        } else {
                            grad[m] = val * alpha;
                        }
                    }
                } else {
                    for (size_t m = 0; m < n_components; m++) {
                        grad[m] = 4;
                    }
                }
                for (size_t m = 0; m < n_components; m++) {
                    current[m] = current[m] + (grad[m]);
                }
            }
            for (size_t m = 0; m < n_components; m++) {
                head_embedding[j + size_head_embedding *m] = current[m];
            }
            epoch_of_next_negative_sample[i] += n_neg_samples * epochs_per_negative_sample[i];   
        }
        alpha = initial_alpha * (1 - static_cast<double>(static_cast<double>(n) / static_cast<double>(n_epochs)));
        alpha4 = alpha * 4;
        alphaNeg4 = alpha * -4;
        if (prArgs!=NULL){
            double nBynTh = floor(fmod((double) n, nTh));
            if (nBynTh == 0) {
                n_epoch = n + 1;
                if (n_epoch < n_epochs) {
                    int continuing=callProgressFcn(prArgs, n_epoch);
                    if (continuing==0){
                        break;
                    }
                }
            }
        }
    }
    mxFree(current);
    mxFree(other);
    mxFree(grad);
    mxFree(sub);
    mxFree(epoch_of_next_negative_sample);
    mxFree(epochs_per_negative_sample);
    mxFree(epoch_of_next_sample);
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    constexpr mwIndex INITIAL_ARGS = 0;
    unsigned verbose = 0;
    if (nlhs != 1){
        if (nrhs<14){
            mexErrMsgTxt("Why call if you don't listen? 0 argout and <14 argin (epochCallbackFcn)");
            return;
        }
    }
/* sanity check: ensure we have input data */
    if (nrhs != 13 && nrhs != 14) {
        mexPrintf("Expecting 13 or 14 args NOT %d args!\n", nrhs);
        mexErrMsgTxt("Expected args: head_embedding, tail_embedding, head, tail, n_epochs, n_vertices, epochs_per_sample, a, b, gamma, initial_alpha, negative_sample_rate, randomize, @(data,epochs)fcnProgress");
        return;
    }

    double *head_embedding, *tail_embedding;
    unsigned *head, *tail;
    double *epochs_per_sample;
    int n_vertices, rand;
    double negative_sample_rate, a, b, gamma, initial_alpha;
    int *randis;
    const mwSize n_components = mxGetN(prhs[INITIAL_ARGS]);
    const mwSize size_head_embedding = mxGetM(prhs[INITIAL_ARGS]);

    mwSize cols=n_components;
    mwSize rows=size_head_embedding;

    plhs[0] = mxDuplicateArray(prhs[INITIAL_ARGS]);
    head_embedding = mxGetPr(plhs[0]);

    mwSize size_tail_embedding = mxGetM(prhs[INITIAL_ARGS + 1]);
    if (verbose>0)
        mexPrintf("Length of tail_embedding is %d!\n", size_tail_embedding);

    mxLogical move_other = size_tail_embedding<1;
    if (verbose>0){
        mexPrintf("move_other=%d && size_tail_embedding=%d\n", move_other, size_tail_embedding);
    }
    int n_epochs = mxGetScalar(prhs[INITIAL_ARGS + 4]);
    int epochReports=10;
    if (!move_other) { //template reduction (supervised or unsupervised)
        mxArray *tailEmbArray = mxDuplicateArray(prhs[INITIAL_ARGS+1]);
        tail_embedding = mxGetPr(tailEmbArray);
        if (size_head_embedding<10001 && n_components<3) {
            epochReports=4;
        } else if (size_head_embedding<20000) {
            epochReports = 5;
        } else if (size_head_embedding<100000){
            if (n_epochs==30)
                epochReports = 6;
            else
                epochReports=10;
        }
    } else {
        tail_embedding = head_embedding;
        size_tail_embedding=size_head_embedding;
        if (size_head_embedding<8194) {
            epochReports = 5;
        } else if (size_head_embedding < 20000){
            epochReports=7;
        }else if (size_head_embedding>90000){
            epochReports = 20;
        }
    }
    const mwSize size_head = mxGetM(prhs[INITIAL_ARGS + 2]);
    if (mxGetClassID(prhs[INITIAL_ARGS + 2]) != mxUINT32_CLASS) {
        mexErrMsgTxt("Oh no! Head array for indexing is wrong data type!");
        return;
    }
    head = static_cast<unsigned *>(mxGetData(prhs[INITIAL_ARGS + 2]));
    double *dHead = mxGetPr(prhs[INITIAL_ARGS + 2]);
    int *testHead = (int *) (dHead);

    const mwSize size_tail =  mxGetM(prhs[INITIAL_ARGS + 3]);
    if (verbose>0)
        mexPrintf("Length of tail is %d!\n", size_tail);

    tail = static_cast<unsigned *>(mxGetData(prhs[INITIAL_ARGS + 3]));

    mxArray *fcnProgress=NULL, *out=NULL, *epochs=NULL;
    int *epochsPtr=NULL;

    if (nrhs==14){
        epochs=mxCreateNumericMatrix(1, 2, mxINT32_CLASS, mxREAL);
        epochsPtr=(int*)mxGetData(epochs);
        epochsPtr[1]=n_epochs;
        fcnProgress=(mxArray *) prhs[13];
        out=(mxArray *) plhs[0];
    } 
    mxArray *args[3]={fcnProgress, out, epochs};
    mxArray **prArgs;
    if (nrhs==14){
        prArgs=args;
    } else {
        prArgs=NULL;
    }
    n_vertices = mxGetScalar(prhs[INITIAL_ARGS + 5]);
    if (verbose>0)
        mexPrintf("n_vertices is %d \n", n_vertices);
    if (move_other && size_head_embedding != n_vertices) {
        mexErrMsgTxt("Whoa! bead_embedding doesn't have n_vertices rows??");
        return;
    }
    if (!move_other && size_tail_embedding != n_vertices) {
        mexErrMsgTxt("Whoa! tail_embedding doesn't have n_vertices rows??");
        return;
    }
    mwSize size_epochs_per_sample = mxGetM(prhs[INITIAL_ARGS + 6]);
    epochs_per_sample = mxGetPr(prhs[INITIAL_ARGS + 6]);

    a = mxGetScalar(prhs[INITIAL_ARGS + 7]);
    if (verbose>0)
        mexPrintf("a value is  %lf \n", a);

    b = mxGetScalar(prhs[INITIAL_ARGS + 8]);

    if (verbose>0)
        mexPrintf("b values is  %lf \n", b);

    gamma = mxGetScalar(prhs[INITIAL_ARGS + 9]);
    if (verbose>0)
        mexPrintf("gamma value is %lf \n", gamma);

    initial_alpha = mxGetScalar(prhs[INITIAL_ARGS + 10]);
    if (verbose>0)
        mexPrintf("initial_alpha values is %lf \n", initial_alpha);

    negative_sample_rate = mxGetScalar(prhs[INITIAL_ARGS + 11]);
    if (verbose>0)
        mexPrintf("negative_sample_rate value is %d \n", negative_sample_rate);
    
    rand = mxGetScalar(prhs[INITIAL_ARGS + 12]);
    if (verbose>0)
        mexPrintf("rand value is %d \n", rand);
    
    if (rand==0)
        srand(0);

    const int n_1_simplices = size_epochs_per_sample;
    if (move_other){
        doMoveOther(head_embedding, tail_embedding, head, tail, n_epochs,
                n_vertices, epochs_per_sample, a, b, gamma, initial_alpha,
                negative_sample_rate, n_1_simplices, n_components,
                prArgs, epochReports);
    } else {
        doNotMoveOther(head_embedding, tail_embedding, head, tail, 
                n_epochs, n_vertices, epochs_per_sample, a, b, gamma, 
                initial_alpha, negative_sample_rate, n_1_simplices,
                n_components, size_head_embedding, prArgs, epochReports);
    }
}

