#ifndef _MEX_ARG_H_
#include <mex.h>
    class mexarg
    {
     public:
         
        int isCellmatrix(int nrhs, const mxArray *prhs[]);
        int require(int nrhs, const mxArray *prhs[], const char *name);
        int getstring(int nrhs, const mxArray *prhs[], const char *name,char *target,int len);
        int getmatrix(int nrhs,const mxArray *prhs[], const char *name,mxArray **target);
        int getscalar(int nrhs,const mxArray *prhs[], const char *name,double *target);
    };
	
#endif