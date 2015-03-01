//--------------------------------------------------------------
// file: fspike.cpp - Matlab Mex function application
//--------------------------------------------------------------
#include <iostream>
#include <math.h>
 
extern "C" {
#include "mex.h"
}

void fspike( int nlhs, mxArray *plhs[], int nrhs, const mxArray  *prhs[] )
{
  
  const double N = mxGetScalar( prhs[0]);
  const double p = mxGetScalar( prhs[1]);
  const double lam = mxGetScalar( prhs[2]);
  const double* const dmu = mxGetPr( prhs[3]);
  const int lent = mxGetM( prhs[3]);
  const int lens = mxGetN( prhs[3]);
  const int lastInd = (lens+1) * (lent+1) - 1;
  
  double* Kn = new double[ lastInd+1];
  double* Kn1 = new double[ lastInd+1];
  for (int i = 0; i < lastInd; i++){
    Kn[ i] = 1;
    Kn1[ i] = 0;
  }
  //  mexPrintf("%f %f %f %i %i %i\n", N, lam, p, lens, lent, lastInd);
  
/* MATLAB-CODE
res = 0;
lam2=lam^2;
for n=1:N
    for s=n:lens
        Kt = 0;
        for t=n:lent
           
            Kt1 = lam2*Kn(s,t)*dmu(s,t) + lam*Kt;
            Kn1(s+1,t+1) = lam*Kn1(s,t+1) + Kt1;
            Kt = Kt1;
        end          
    end
    res = res + Kn1(lens+1,lent+1)*p^n;
    Kn = Kn1;
    Kn1 = zeros(lens+1, lent+1);
end
*/
 
  double res = 0;
  double Kt,Kt1;
  const double lam2 = lam *lam;
  double pp = p;
  int ss, ss1;

  for (int n = 0; n<N; n++){
    for (int s = n; s<lens; s++){
      Kt = 0;
      ss = s*(lent+1);
      ss1 = ss + lent + 1;
      for (int t = n; t<lent; t++){
	
	Kt1 = lam2 * Kn[ ss + t] * dmu[ s*lent + t] + lam*Kt;
	Kn1[ ss1 + t+1] = lam * Kn1[ ss + t+1] + Kt1;
	//    mexPrintf( "%e %e\n", Kt1, Kn1[ ss1 + t+1]);
	//    mexPrintf( "%e \n", Kn[ ss+t]);
	Kt = Kt1;
      }
    }
    res = res + Kn1[ lastInd] *pp;
    pp = pp*p;
    for (int i = 0; i < lastInd; i++){
      Kn[ i] = Kn1[ i];
      Kn1[ i] = 0;
    }
  }
  
  delete[] Kn;
  delete[] Kn1;

  plhs[0]= mxCreateDoubleMatrix(1, 1, mxREAL);
  double * result=mxGetPr(plhs[0]);
  *result=res;
  
  
} // end fspike()

extern "C" {
  //--------------------------------------------------------------
  // mexFunction - Entry point from Matlab. From this C function,
  //   simply call the C++ application function, above.
  //--------------------------------------------------------------
  void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray  *prhs[] )
  {
    fspike(nlhs, plhs, nrhs, prhs);
  }
}
