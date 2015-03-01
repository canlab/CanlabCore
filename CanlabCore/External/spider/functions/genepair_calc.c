

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void mexFunction (int nlhs, mxArray *plhs[], int nrhs,
                  const mxArray *prhs[])
{
  /* mxGetString (prhs [0], fn, mxGetN(prhs [0]) + 1); 
     mexEvalString("name=[];");*/

  double *X,*Y,*K;
  int l,n,l2,n2;  
  int p,q,i; double xi,yi,xi1,yi1,k,f1,f2,f3;

  l = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  X = (double *) mxGetPr(prhs[0]);
  l2 = mxGetM(prhs[1]);
  n2 = mxGetN(prhs[1]);
  Y = (double *) mxGetPr(prhs[1]);

  /*  printf("i am calculating the kernel on %dx%d and %dx%d...\n",l,n,l2,n2); */

  plhs[0] = mxCreateDoubleMatrix(l,l2,mxREAL);  
  K=mxGetPr(plhs[0]);
  
  /*
  for(p=0;p<l;p++)
 {
   for(q=0;q<n;q++)
     {
          printf("%f ",X[q*l+p]);
     }
   printf("\n");
   }*/

  for(p=0;p<l;p++)
  for(q=0;q<l2;q++)
  {
    k=0; f1=0; f2=0; f3=0;
    for(i=0;i<n;i++)
      {
	xi=X[p+l*i];    yi=Y[q+l2*i];
	xi1=X[p+l*(i-1)]; yi1=Y[q+l2*(i-1)];
        k=k+(xi*yi)*(i);
        if (i>0) 
	  {
	    f1=f1+yi1;    
            f2=f2+xi1;    
            f3=f3+xi1*yi1;
	  }
        k=k-f1*xi-yi*f2+f3; 
      }
    K[p+l*q]=k;
  }


  /*  for(k=0;k<4;k++)
    R[i+k*rows]=num[k];
    i++;*/

  nlhs=1; /* return one thing, (the matrix) */
}






