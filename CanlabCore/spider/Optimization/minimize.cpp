#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mex.h>
#include "VectorUtil.hpp"

#ifndef min
template <class T> inline T min(T x,T y) { return (x<y)?x:y; }
#endif
#ifndef max
template <class T> inline T max(T x,T y) { return (x>y)?x:y; }
#endif


/* % Minimize a continuous differentialble multivariate function. Starting point
% is given by "X" (D by 1), and the function named in the string "f", must
% return a function value and a vector of partial derivatives. The Polack-
% Ribiere flavour of conjugate gradients is used to compute search directions,
% and a line search using quadratic and cubic polynomial approximations and the
% Wolfe-Powell stopping criteria is used together with the slope ratio method
% for guessing initial step sizes. Additionally a bunch of checks are made to
% make sure that exploration is taking place and that extrapolation will not
% be unboundedly large. The "length" gives the length of the run: if it is
% positive, it gives the maximum number of line searches, if negative its
% absolute gives the maximum allowed number of function evaluations. You can
% (optionally) give "length" a second component, which will indicate the
% reduction in function value to be expected in the first line-search (defaults
% to 1.0). The function returns when either its length is up, or if no further
% progress can be made (ie, we are at a minimum, or so close that due to
% numerical problems, we cannot get any closer). If the function terminates
% within a few iterations, it could be an indication that the function value
% and derivatives are not consistent (ie, there may be a bug in the
% implementation of your "f" function). The function returns the found
% solution "X", a vector of function values "fX" indicating the progress made
% and "i" the number of iterations (line searches or function evaluations,
% depending on the sign of "length") used.
%
% Usage: X = minimize(X, f, length, P)
%
% See also: checkgrad 
%
% Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13*/

double* evaluate(const char* fname, double* X, int nX, const mxArray* P, double* fX){	
	mxArray* lhs[2];
	mxArray* rhs[2];
	rhs[0] = mxCreateNumericMatrix(nX, 1, mxDOUBLE_CLASS, mxREAL);	
	memcpy(mxGetPr(rhs[0]), X, nX * mxGetElementSize(rhs[0]));
	rhs[1] = (mxArray*)P;	
	mexCallMATLAB(2, lhs, 2, rhs, fname);	// [fX grad]=f(X,P)	
	*fX = mxGetScalar(lhs[0]);
	double* grad = mxGetPr(lhs[1]);	
	mxDestroyArray(rhs[0]);		
	mxDestroyArray(lhs[0]);		
// 	mexPrintf(".");
// 	fflush(stdout);
	return grad;
}


double* minimize(double* X, int nX, const char* fname, int length, const mxArray* P, double* fX){	
	const double RHO=0.01;                 //           % a bunch of constants for line searches
	const double SIG=0.5;     //  % RHO and SIG are the constants in the Wolfe-Powell conditions
	const double INT=0.1;    //% don't reevaluate within 0.1 of the limit of the current bracket
	const double EXT=3.0;      //              % extrapolate maximum 3 times the current bracket/*
	const int MAXEVAL=20;         //                % max 20 function evaluations per line search
	const double RATIO=100.0;        //                              % maximum allowed slope ratio
	const double MINDOUBLE=2.2251e-308;
	
	char* S1="Linesearch";
	char* S2="Function evaluation"; 
	char* S = (length > 0)? S1:S2;	
	
	double* X0;
    	double f0;
    	double* df0;
    	double red = 1;
    	int iterations_i = 0;  
    	int M;
    	bool ls_failed = false;
    	double f1,f2;
    	double* df2;
    	double d2, f3, d3, z3, A, B, limit, z2;
    	bool success;
		    
	/*double f1 = f.getValue(X); 
	double[] df1 = f.getGradient(X); // partial dervatives*/
	double* df1 = evaluate(fname, X, nX, P, &f1);	
	iterations_i = iterations_i + (length<0);                 //count epochs?!	
	double* s = svmult(-1.0, df1, nX); //search direction is steepest	
	double d1 = 0;
	for (int i = 0; i < nX; i++) 
		d1 += -s[i] * s[i]; //this is the slope
	double z1 = red / (1.0 - d1); // initial step is red/(|s|+1)
	while (iterations_i  < abs(length)) {	    		      	      
	      iterations_i = iterations_i + (length>0);  //count iterations?!
	      X0 = clone(X,nX);
	      f0 = f1;
	      df0 = clone(df1, nX);//make a copy of current values
	      X = vvadd(X, svmult(z1, s, nX), nX);        // begin line search      
	      
	      /*f2 = f.getValue(X);      
	      df2 = f.getGradient(X); //and gradient*/
	      df2 = evaluate(fname, X, nX, P, &f2);
	      iterations_i = iterations_i + (length<0);                 //count epochs?!	      
	      
	      d2 = scalar(df2, s, nX); 
	      f3 = f1; 
	      d3 = d1; 
	      z3 = -z1;             // initialize point 3 equal to point 1	      	      
	      M = (length > 0)? MAXEVAL : (min(MAXEVAL, -length - iterations_i));	      
	      	      
	      success = false;
	      limit = -1.0;                    // initialize quanteties
	      z2 = 0;
	      while (true) {
		        while (((f2 > f1+z1*RHO*d1) || (d2 > -SIG*d1)) && (M > 0)) {
			          limit = z1;
			          if (f2 > f1) {
			            z2 = z3 - (0.5*d3*z3*z3)/(d3*z3+f2-f3);                 // quadratic fit
			          } 
			          else {
			            A = 6.0*(f2-f3)/z3+3.0*(d2+d3);                                // cubic fit
			            B = 3.0*(f3-f2)-z3*(d3+2.0*d2);
			            z2 = (sqrt(B*B-A*d2*z3*z3)-B)/A;      // numerical error possible - ok!
			          }
			          if (isinf(z2) || isnan(z2)) 
			          	z2 = z3/2.0;
			          z2 = max(min(z2, INT*z3),(1.0-INT)*z3);  // don't accept too close to limits
			          z1 = z1 + z2;                                           // update the step
			          X = vvadd(X, svmult(z2,s,nX),nX);          			          	          	  
			          
				  /*f2 = f.getValue(X);		         
			          df2 = f.getGradient(X);*/
				  df2 = evaluate(fname, X, nX, P, &f2);
				  M = M -1;
				  iterations_i = iterations_i + (length<0);                 //count epochs?!
				  
			          d2 = scalar(df2, s, nX);
			          z3 = z3-z2;          
		        } // end of while
		        if ((f2 > f1+z1*RHO*d1) || (d2 > -SIG*d1)) {
		        	break;                           // this is a failure
			} else if (d2 > SIG*d1) {
				success = true;
				break;                                             // success
		        }
		        else if (M == 0) {
		        	break;
		        }
		        A = 6.0*(f2-f3)/z3+3.0*(d2+d3);                      // make cubic extrapolation
		        B = 3.0*(f3-f2)-z3*(d3+2*d2);
		        z2 = -d2*z3*z3/(B+sqrt(B*B-A*d2*z3*z3));        // num. error possible - ok!
		        if  (isnan(z2) || isinf(z2) || (z2 < 0)) {   // num prob or wrong sign?
				if (limit < -0.5) {                             // if we have no upper limit
					z2 = z1 * (EXT-1);                           // the extrapolate the maximum amount
				} else {
					z2 = (limit-z1)/2;                           // otherwise bisect
				}
		        } 
		        else if ((limit > -0.5) & (z2+z1 > limit)) {     // extraplation beyond max?
		        	z2 = (limit-z1)/2;    // bisect
		        } else if ((limit < -0.5) & (z2+z1 > z1*EXT)) {       // extrapolation beyond limit
		            	z2 = z1*(EXT-1.0);                           // set to extrapolation limit
		        } else if (z2 < -z3*INT) {
		            	z2 = -z3*INT;
		        } else if ((limit > -0.5) & (z2 < (limit-z1)*(1.0-INT)) )  { // too close to limit?
		            	z2 = (limit-z1)*(1.0-INT);
		        }
		
		        f3 = f2;
		        d3 = d2;
		        z3 = -z2;                  // set point 3 equal to point 2
		        z1 = z1 + z2;
		        X = vvadd(X, svmult(z2,s,nX),nX);                     // update current estimates	
			/*	
		        f2 = f.getValue(X);
		        df2 = f.getGradient(X);
			*/
			df2 = evaluate(fname, X, nX, P, &f2);
			M = M -1; 
			iterations_i = iterations_i + (length<0);                 //count epochs?!
			
		        d2 = scalar(df2, s, nX);		                
	        }     
	
	
		if (success) {                                        // if line search succeeded		            
			f1 = f2;
			mexPrintf("%s %6i;  Value %4.6e\r", S, iterations_i, f1);
			fflush(stdout);
			*fX = f1;
			s = vvadd(svmult((scalar(df2,df2,nX) - scalar(df1,df2,nX)) / (scalar(df1,df1,nX)),s,nX), svmult(-1.0,df2,nX),nX);     // Polack-Ribiere direction   
			double* tmp = df1;
			df1 = df2;
			df2 = tmp;                         // swap derivatives
			d2 = scalar(df1,s,nX);
			if (d2 > 0) {
				s = svmult(-1.0, df1, nX);                              // otherwise use steepest direction
				d2 = scalar(svmult(-1.0,s,nX), s, nX);
			}
			z1 = z1 * min(RATIO, d1/(d2-MINDOUBLE));          // slope ratio but max RATIO
			d1 = d2;
			ls_failed = false;                              // this line search did not fail   
		} 
		else {
			X = X0;
			f1 = f0;
			df1 = df0;
			if (ls_failed || (iterations_i > abs(length))) 
				break;        // line search failed twice in a row, or we ran out of time, so we give up
			double* tmp = df1;
			df1 = df2;
			df2 = tmp;                         // swap derivatives
			s = svmult(-1.0, df1, nX);                          // try steepest
			d1 = scalar(svmult(-1.0, s, nX), s, nX);
			z1 = 1.0/(1.0-d1);	   	    
			ls_failed = true;                                    // this line search failed
		}   	    		
 	}
	mexPrintf("\n");
	return X;
}

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]){
        
    if((nrhs < 4) || (nlhs > 2)){
    	mexErrMsgTxt("Usage: [X fX]=minimize(X,f,length,P)");
        return;
    }    
    double* X = mxGetPr(prhs[0]);
    int m = mxGetM(prhs[0]);
    int n = mxGetN(prhs[0]);
    int nX = max(m,n);
    int flen = mxGetN(prhs[1])+1;    
    char* fname = (char*) mxCalloc(flen, sizeof(char));
    mxGetString(prhs[1], fname, flen);    
    int length = (int)mxGetScalar(prhs[2]);
    const mxArray* P;
    if(nrhs == 4)
    	P = prhs[3];
    double fX;
    X = minimize(X, nX, fname, length, P, &fX);
    plhs[0] = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);	
    memcpy(mxGetPr(plhs[0]), X, nX * mxGetElementSize(plhs[0]));    
    if (nlhs == 2)
    	plhs[1] = mxCreateDoubleScalar(fX);
}
