#include <mex.h>
#include <stdlib.h>
#include <string.h>
#include <mexarg.h>

extern "C"{
	
    
    
	void mexFunction(
		int nlhs,              // Number of left hand side (output) arguments
		mxArray *plhs[],       // Array of left hand side arguments
		int nrhs,              // Number of right hand side (input) arguments
		const mxArray *prhs[]  // Array of right hand side arguments
		)
	{
        mexarg arg;
        
		double c=0;
        mxArray *x;
        int *dims;
        int ndims;
        char string[128]="";
      
        int result=1;
        
        
        
        if (!(arg.require(nrhs,prhs,"C") & 
            arg.require(nrhs,prhs,"x")))
        {
            mexPrintf("We require parameters C and x!\n");
        }
        else
        {
            
            result&=arg.getscalar(nrhs,prhs,"C",&c);
            result&=arg.getmatrix(nrhs,prhs,"x",&x);
            result&=arg.getstring(nrhs,prhs,"n",string,128);
        
            mexPrintf("%s\n",string);
            mexPrintf("%d\n",mxGetM(x));
            mexPrintf("%f\n",c);
            
           }
	}
	
}
