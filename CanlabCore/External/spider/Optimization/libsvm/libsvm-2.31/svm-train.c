#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "svm.h"
#include "mex.h"
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))



struct svm_model
{
	struct svm_parameter param;	// parameter
	int nr_class;		// number of classes, = 2 in regression/one class svm
	int l;			// total #SV
	struct svm_node **SV;		// SVs (SV[l])
	double **sv_coef;	// coefficients for SVs in decision functions (sv_coef[n-1][l])
	double *rho;		// constants in decision functions (rho[n*(n-1)/2])

	// for classification only

	int *label;		// label of each class (label[n])
	int *nSV;		// number of SVs for each class (nSV[n])
				// nSV[0] + nSV[1] + ... + nSV[n-1] = l
	// XXX
	int free_sv;		// 1 if svm_model is created by svm_load_model
				// 0 if svm_model is created by svm_train
};



void exit_with_help()
{
	printf(
	"Usage: svm-train [options] training_set_file [model_file]\n"
	"options:\n"
	"-s svm_type : set type of SVM (default 0)\n"
	"	0 -- C-SVC\n"
	"	1 -- nu-SVC\n"
	"	2 -- one-class SVM\n"
	"	3 -- epsilon-SVR\n"
	"	4 -- nu-SVR\n"
	"-t kernel_type : set type of kernel function (default 2)\n"
	"	0 -- linear: u'*v\n"
	"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
	"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	"-d degree : set degree in kernel function (default 3)\n"
	"-g gamma : set gamma in kernel function (default 1/k)\n"
	"-r coef0 : set coef0 in kernel function (default 0)\n"
	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
	"-m cachesize : set cache memory size in MB (default 40)\n"
	"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
	"-h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
	"-wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
	"-v n: n-fold cross validation mode\n"
	);
	exit(1);
}

void do_cross_validation();

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
int cross_validation = 0;
int nr_fold;

void do_cross_validation()
{
	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;

	// random shuffle
	for(i=0;i<prob.l;i++)
	{
		int j = rand()%(prob.l-i);
		struct svm_node *tx;
		double ty;
			
		tx = prob.x[i];
		prob.x[i] = prob.x[j];
		prob.x[j] = tx;

		ty = prob.y[i];
		prob.y[i] = prob.y[j];
		prob.y[j] = ty;
	}

	for(i=0;i<nr_fold;i++)
	{
		int begin = i*prob.l/nr_fold;
		int end = (i+1)*prob.l/nr_fold;
		int j,k;
		struct svm_problem subprob;

		subprob.l = prob.l-(end-begin);
		subprob.x = Malloc(struct svm_node*,subprob.l);
		subprob.y = Malloc(double,subprob.l);
			
		k=0;
		for(j=0;j<begin;j++)
		{
			subprob.x[k] = prob.x[j];
			subprob.y[k] = prob.y[j];
			++k;
		}
		for(j=end;j<prob.l;j++)
		{
			subprob.x[k] = prob.x[j];
			subprob.y[k] = prob.y[j];
			++k;
		}

		if(param.svm_type == EPSILON_SVR ||
		   param.svm_type == NU_SVR)
		{
			struct svm_model *submodel = svm_train(&subprob,&param);
			double error = 0;
			for(j=begin;j<end;j++)
			{
				double v = svm_predict(submodel,prob.x[j]);
				double y = prob.y[j];
				error += (v-y)*(v-y);
				sumv += v;
				sumy += y;
				sumvv += v*v;
				sumyy += y*y;
				sumvy += v*y;
			}
			svm_destroy_model(submodel);
			printf("Mean squared error = %g\n", error/(end-begin));
			total_error += error;			
		}
		else
		{
			struct svm_model *submodel = svm_train(&subprob,&param);
			int correct = 0;
			for(j=begin;j<end;j++)
			{
				double v = svm_predict(submodel,prob.x[j]);
				if(v == prob.y[j])
					++correct;
			}
			svm_destroy_model(submodel);
			printf("Accuracy = %g%% (%d/%d)\n", 100.0*correct/(end-begin),correct,(end-begin));
			total_correct += correct;
		}

		free(subprob.x);
		free(subprob.y);
	}		
	if(param.svm_type == EPSILON_SVR || param.svm_type == NU_SVR)
	{
		printf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
		printf("Cross Validation Squared correlation coefficient = %g\n",
			((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
			((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
			);
	}
	else
		printf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/prob.l);
}



//////////////////////////////////////////////////////////

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    
    int elements, max_index, i, j;
    int m=0,d=0;
    int ifield;
    unsigned int neq=0;
    double valmax=1e50;
    double *xsv, *alpha, *Xtmp, b, epsilon, C;
    int nb_rows, n;

    if (nrhs > 3 || nrhs < 2) {
	    mexErrMsgTxt("Usage: [x,how] "
		         "= svm-train(X,Y,options)");
	    return;
    }
    switch (nrhs) {
    case 2:
	     {
		    if ((!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		     ||  mxIsSparse(prhs[0]))||(!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
		     ||  mxIsSparse(prhs[1])))
		      {
			 mexErrMsgTxt("Datasets must be "
				      "a real.");
			 return;
		    }
		   }
    case 3:
	    if (mxGetM(prhs[5]) != 0 || mxGetN(prhs[5]) != 0) {
		    if (!mxIsStruct(prhs[2]))
		    {
			 mexErrMsgTxt("options must be "
				      "a struct.");
			 return;
		    }
	
	    }
	    }
	    
	    //////////////////////////
    
    if (nlhs > 5 || nlhs < 1) {
	    mexErrMsgTxt("Usage: [alpha,b,Xsv,epsilon,C] "
		         "= train-svm(X,Y,options)");
	    return;
    }

    switch (nlhs) {
    case 5:
    case 4:
    case 3:
        m = mxGetM(prhs[0]);
        d = mxGetN(prhs[0]);
        plhs[2]=mxCreateDoubleMatrix(m,d,mxREAL); // for Xsv
        xsv = mxGetPr(plhs[2]);
    case 2:
    case 1:
	    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL); // for alpha
	    alpha = mxGetPr(plhs[0]);
	    break;
    }
	
    
    // Read the values from matlab to C
    
    // for the dataset
	prob.l = m;
	elements = m*d;
	Xtmp = mxGetPr(prhs[0]);
	
	prob.y = Malloc(double,m);
	prob.x = Malloc(struct svm_node *,m);
	x_space = Malloc(struct svm_node,elements);

	max_index = d;
	j=0;
	prob.y = mxGetPr(prhs[1]);
	for(i=0;i<m;i++)
	{
	  for(j=0;j<d;j++)
	  {	  
		prob.x[i] = &x_space[j];
		x_space[j].index=j;
		x_space[j].value=Xtmp[i*d+j];
		}	
	}

	// for param
	// default values
	
	
	param.svm_type = mxGetScalar(mxGetField(prhs[2],1,"svm_type"));
	param.kernel_type = mxGetScalar(mxGetField(prhs[2],1,"kernel_type"));
	param.degree = mxGetScalar(mxGetField(prhs[2],1,"degree"));
	param.gamma = mxGetScalar(mxGetFieldByNumber(prhs[2],0,4));	// 1/k
	param.coef0 = mxGetScalar(mxGetFieldByNumber(prhs[2],0,5));
	param.nu = mxGetScalar(mxGetFieldByNumber(prhs[2],0,6));
	param.cache_size = mxGetScalar(mxGetFieldByNumber(prhs[2],0,7));
	param.C = mxGetScalar(mxGetFieldByNumber(prhs[2],0,8));
	param.eps = mxGetScalar(mxGetFieldByNumber(prhs[2],0,9));
	param.p = mxGetScalar(mxGetFieldByNumber(prhs[2],0,10));
	param.shrinking = mxGetScalar(mxGetFieldByNumber(prhs[2],0,11));
	param.nr_weight = mxGetScalar(mxGetFieldByNumber(prhs[2],0,12));
	param.weight_label =NULL;// mxGetScalar(mxGetFieldByNumber(prhs[2],0,13));
	param.weight = NULL;//mxGetScalar(mxGetFieldByNumber(prhs[2],0,14));
	
	
	// put options directly to param
	
	if(param.gamma == 0)
		param.gamma = 1.0/max_index;

		
		
    if(cross_validation)
	{
		do_cross_validation();
	}
	else
	{
		model = svm_train(&prob,&param);
	}

	
	/// build the output for matlab
	
	
	free(prob.y);
	free(prob.x);
	free(x_space);
	
	/// from the model to...
	for(i=0;i<model->l;i++){
	alpha[i] = model->sv_coef[1][i];
	xsv[i] = model->SV[i][j].value;
	}
	plhs[1] = model->rho;
	plhs[3] = &model->param.p;
	plhs[4] = &model->param.C;
}

