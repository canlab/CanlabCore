//#include "svqp.h"

#include <stdlib.h>
#include <math.h>


// svreal - typedef for floating point type (float or double)
typedef double svreal;

// svbool - typedef for boolean value
typedef char svbool;


/* n -- dimension of vectors x,cmin,cmax */

int n;

// cmin, cmax -- problem parameters
  svreal *cmin, *cmax;
  // valmax -- maximum value of the objective function
  svreal valmax;
  // sumflag -- set if linear constraint must be considered
  svbool sumflag;     
  // ktflag -- set if optimizer must loop until all kt conditions are met.
  svbool ktflag;     
  // x -- current solution vector
  svreal *x;
  // w -- current solution value
  svreal w;
  // error -- description of error
  char *error;
  // epsgr -- stopping criterion on gradient norm
  svreal epsgr;
  // epskt -- clamping proximity for constraints
  svreal epskt;
  // maxst -- maximal gradient step value
  svreal maxst;

  int        iterations;
  svreal    *g;
  svreal    *z;
  svreal     gnorm;
  // info about past iteration
  svbool     restartp;
  svreal    *gsav;
  svreal     zgsav;
  // clamping information
  char      *clamp;
  int        nactive;
  // memory allocation
  svreal    *mem;
  svreal *tmp;
  // b -- linear term
  svreal *b;
  // a -- quadratic matrix (line order)
  svreal *a;

// ================================================================
// Utilities
// ================================================================

svreal dot(svreal *x, svreal *y)
{
  int i;
  svreal sum=0.0;
  for (i=0; i<n; i++)
    sum += *x++ * *y++;
  return sum;
}


#define max(a,b) ((a) > (b)) ? (a) : (b)

// ================================================================
// SimpleQuadraticProgram
// ================================================================
//     Maximize          - x.A.x + b.x
//     with              Cmin_i <= x_i <= Cmax_i
//     and (optionally)  sum_i x_i = 0


void init()
{
  // allocate memory
  mem     = (svreal *)malloc(6*n*sizeof(svreal));
  tmp     = (svreal *)malloc(n*sizeof(svreal));
  clamp   = (char *)malloc(n*sizeof(char));
  x       = mem;
  cmin    = mem+n;
  cmax    = mem+2*n;
  g       = mem+3*n;
  z       = mem+4*n;
  gsav    = mem+5*n;
  // initialize variables
  valmax  = 1e20;
  sumflag = 0;
  ktflag  = 1;
  w       = 0;
  error   = 0;
  epsgr   = (svreal)1E-8;
  epskt   = (svreal)1E-20;
  maxst   = (svreal)1E+10;
  nactive = n;
}

void free_mem()
{
  free(mem);
  free(clamp);
}

void compute_Ax(svreal *x, svreal *y)
{
  int i;
  svreal *aa = a;
  for (i=0; i<n; i+=1,aa+=n)
    y[i] = dot(x,aa);
}

svreal compute_gx(svreal *x, svreal *g)
{
  int i;
  compute_Ax(x,g);
  for (i=0; i<n; i++)
    {
      tmp[i] = g[i]/2 - b[i];
      g[i] = g[i] - b[i];
    }
  return dot(tmp,x);
}

svreal compute_ggx(svreal *z)
{
  compute_Ax(z,tmp);    // can be faster when A is symetric
  return dot(tmp,z);
}

// ================================================================
// ConvexProgram
// ================================================================


#define CLAMPNO  0
#define CLAMPMIN 1
#define CLAMPMAX 2

void project_with_linear_constraint(void)
{
  int i;
  int oldactive;
  svreal gcorr;
  // Iterate (maximum n times, usually two times)
  do 
  {
    oldactive = nactive;
    // Compute gradient correction
    gcorr = 0;
    for (i=0; i<n; i++)
      {
        if (clamp[i] == CLAMPNO)
          gcorr += g[i];
        else
          g[i] = 0;
      }
    gcorr = gcorr / nactive;
    // Check constraints
    for (i=0; i<n; i++)
      {
        if (clamp[i] == CLAMPNO)
          {
            g[i] = g[i] - gcorr;
            if (g[i] < 0)
              {
                if (x[i] <= cmin[i]+epskt)
                  {
                    x[i]      = cmin[i];
                    g[i]      = 0;
                    clamp[i]  = CLAMPMIN;
                    restartp  = 1;
                    nactive   -= 1;
                  }
              }
            else
              {
                if (x[i] >= cmax[i]-epskt)
                  {
                    x[i]      = cmax[i];
                    g[i]      = 0;
                    clamp[i]  = CLAMPMAX;
                    restartp  = 1;
                    nactive   -= 1;
                  }
              }
          }
      }
  } 
  while (nactive < oldactive);
}

void project_without_linear_constraint(void)
{
  int i;
  for (i=0; i<n; i++)
    {
      if (clamp[i] != CLAMPNO)
        {
          g[i] = 0;
        }
      else if (g[i] < 0)
        {
          if (x[i] <= cmin[i]+epskt)
            {
              x[i]      = cmin[i];
              g[i]      = 0;
              clamp[i]  = CLAMPMIN;
              restartp  = 1;
              nactive   -= 1;
            }
        }
      else
        {
          if (x[i] >= cmax[i]-epskt)
            {
              x[i]      = cmax[i];
              g[i]      = 0;
              clamp[i]  = CLAMPMAX;
              restartp  = 1;
              nactive   -= 1;
            }
        }
    }
}

svbool adjust_clamped_variables(void)
{
  int i;
  int oldactive = nactive;
  // Recompute gradient again
  w = - compute_gx(x, g);
  for (i=0; i<n; i++)
    g[i] = - g[i];
  // Reset clamp status
  nactive = n;
  for (i=0; i<n; i++)
    clamp[i] = CLAMPNO;
  // Project gradient and recompute clamp status
  if (sumflag)
    project_with_linear_constraint();
  else
    project_without_linear_constraint();
  // Recompute gradient_norm
  gnorm = dot(g,g);
  // Return true if variables have been unclamped
  return (nactive > oldactive);
}

int run(void)
{
  // initialize variables
  int i;
  int itercg = 0;
  svreal avgw = 0;
  svreal beta;
  svreal step, ostep, curvature;
  for (i=0; i<n; i++)
    clamp[i] = CLAMPNO;
  restartp = 1;
  iterations = 0;
  // test feasibility
  for (i=0; i<n; i++)
    if (x[i]<cmin[i] || x[i]>cmax[i])
      {
        printf("initial x vector is not feasible : x %f  cmin %f  cmax %f\n",x[i], cmin[i], cmax[i]);
        return -1;
      }
  // averaged w for detecting instabilities

  // main loop
  for (;;)
    {
      printf("new iteration : %d\n", iterations);
      // compute gradient
      w = - compute_gx(x, g);
      if (w > valmax) {
        printf("Objective function exceeds valmax");
        return -1;
      }
      for (i=0; i<n; i++) {
        g[i] = - g[i];
      }
      // test instabilities
      if (iterations < n)
        avgw = w;
      else
        avgw += (w-avgw)/(5*n);
      if (w<avgw && iterations>n+n)
        {
          printf("Numerical instability (EPSGR is too small)");
          return -1;
        }
      // project gradient and compute clamp status
      if (sumflag)
        project_with_linear_constraint();
      else
        project_without_linear_constraint();
      // compute gradient_norm and test termination
      gnorm = dot(g,g);
      if (gnorm<epsgr)
        {
          if (! ktflag)
            break;
          // Reexamine clamping status
          if (! adjust_clamped_variables())
            break;
          if (gnorm<epsgr)
            break;
          // Continue processing
        }
      // compute search direction
      if (restartp)
        {
          // Just copy gradient into search direction
          itercg = 0;
          for (i=0; i<n; i++)
            z[i] = gsav[i] = g[i];
          // Save previous gradient
          restartp = 0;
          zgsav = gnorm;
        }
      else
        {
          // Self restarting Hestenes Stiefel
          for (i=0; i<n; i++)
            gsav[i] = g[i] - gsav[i];
          beta = dot(g,gsav) / zgsav;
          for(i=0; i<n; i++)
            z[i] = g[i] + beta * z[i];
          // Save previous gradient
          restartp = 0;
          for (i=0; i<n; i++)
            gsav[i] = g[i];
          zgsav = dot(z,gsav);
        }
      // Compute clipping constraint
      step = maxst;
      for (i=0; i<n; i++)
        {
          if (z[i] < 0)
            {
              ostep = (cmin[i]-x[i]) / z[i];
              if (ostep < step)
                step = ostep;
            }
          else if (z[i] > 0)
            {
              svreal ostep = (cmax[i]-x[i]) / z[i];
              if (ostep < step)
                step = ostep;
            }
        }
      // Compute optimal quadratic step
      curvature = compute_ggx(z);
      if (curvature >= epskt)
        {
          ostep = zgsav / curvature;
          if (ostep < step)
            step = ostep;
        }
      else if (step >= maxst)
        {
          printf("Function is not convex (negative curvature)");
          return -1;
        }
      // perform update
      for (i=0; i<n; i++)
        x[i] += step * z[i];
      // count iterations
      iterations += 1;
      itercg += 1;
    }
  // Termination
  return iterations;
}

//int main(){
//}

void svqp(int *nn, double *Q, double *c, double *l, double *u, double *x0, int *neq, double *vm, double *vect, int *res) {

  int i,j;

  n = *nn; 
  init();
  a = Q;
  b = c;
  for (i=0; i<n; i++) { 
    x[i] = x0[i];
    cmin[i] = l[i];
    cmax[i] = u[i];
  }

  sumflag = *neq;
  valmax = *vm;

  *res = run();

  for (i=0; i<n; i++) { 
    vect[i] = x[i];
  }

  free_mem();
}

void mexFunction(
    int nlhs,  mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
)
{
    double *c=NULL, *Q=NULL, 
	   *l=NULL, *u=NULL, *x=NULL, *x0=NULL;
    int n=0;
    unsigned int neq=0;
    double valmax=1e50;

    if (nrhs > 7 || nrhs < 5) {
	    mexErrMsgTxt("Usage: [x,how] "
		         "= sqvp(H,c,l,u,x0,eq,vm)");
	    return;
    }
    switch (nrhs) {
    case 7:
	    if (mxGetM(prhs[6]) != 0 || mxGetN(prhs[6]) != 0) {
		    if (!mxIsNumeric(prhs[6]) || mxIsComplex(prhs[6]) 
		     ||  mxIsSparse(prhs[6])
		     || !(mxGetM(prhs[6])==1 && mxGetN(prhs[6])==1)) {
			 mexErrMsgTxt("Seventh argument (valeur max) must be "
				      "a real.");
			 return;
		    }
		    valmax = *mxGetPr(prhs[6]);
	    }
    case 6:
	    if (mxGetM(prhs[5]) != 0 || mxGetN(prhs[5]) != 0) {
		    if (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5]) 
		     ||  mxIsSparse(prhs[5])
		     || !(mxGetM(prhs[5])==1 && mxGetN(prhs[5])==1)) {
			 mexErrMsgTxt("Sixth argument (eq) must be "
				      "a boolean.");
			 return;
		    }
		    neq = (unsigned int) *mxGetPr(prhs[5]);
	    }
    case 5:
	    if (mxGetM(prhs[4]) != 0 || mxGetN(prhs[4]) != 0) {
		    if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) 
		     ||  mxIsSparse(prhs[4])
		     || !mxIsDouble(prhs[4]) 
		     ||  mxGetN(prhs[4])!=1 ) {
			 mexErrMsgTxt("Fifth argument (x0) must be "
				      "a column vector.");
			 return;
		    }
		    x0 = mxGetPr(prhs[4]);
		    n  = mxGetM(prhs[4]);
	    }
    case 4:
	    if (mxGetM(prhs[3]) != 0 || mxGetN(prhs[3]) != 0) {
		    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) 
		     ||  mxIsSparse(prhs[3])
		     || !mxIsDouble(prhs[3]) 
		     ||  mxGetN(prhs[3])!=1 ) {
			 mexErrMsgTxt("Fourth argument (u) must be "
				      "a column vector.");
			 return;
		    }
		    if (n != 0 && n != mxGetM(prhs[3])) {
			 mexErrMsgTxt("Dimension error (arg 4 and later).");
			 return;
		    }
		    u = mxGetPr(prhs[3]);
		    n = mxGetM(prhs[3]);
	    }
    case 3:
	    if (mxGetM(prhs[2]) != 0 || mxGetN(prhs[2]) != 0) {
		    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) 
		     ||  mxIsSparse(prhs[2])
		     || !mxIsDouble(prhs[2]) 
		     ||  mxGetN(prhs[2])!=1 ) {
			 mexErrMsgTxt("Third argument (l) must be "
				      "a column vector.");
			 return;
		    }
		    if (n != 0 && n != mxGetM(prhs[2])) {
			 mexErrMsgTxt("Dimension error (arg 3 and later).");
			 return;
		    }
		    l = mxGetPr(prhs[2]);
		    n = mxGetM(prhs[2]);
	    }
    case 2:
	    if (mxGetM(prhs[1]) != 0 || mxGetN(prhs[1]) != 0) {
		    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) 
		     ||  mxIsSparse(prhs[1])
		     || !mxIsDouble(prhs[1]) 
		     ||  mxGetN(prhs[1])!=1 ) {
			 mexErrMsgTxt("Second argument (c) must be "
				      "a column vector.");
			 return;
		    }
		    if (n != 0 && n != mxGetM(prhs[1])) {
			 mexErrMsgTxt("Dimension error (arg 2 and later).");
			 return;
		    }
		    c = mxGetPr(prhs[1]);
		    n = mxGetM(prhs[1]);
	    }
    case 1:
	    if (mxGetM(prhs[0]) != 0 || mxGetN(prhs[0]) != 0) {
		    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) 
		     ||  mxIsSparse(prhs[0]) ) {
			 mexErrMsgTxt("First argument (H) must be "
				      "a matrix.");
			 return;
		    }
		    if (n != 0 && n != mxGetM(prhs[0])) {
			 mexErrMsgTxt("Dimension error (arg 1 and later).");
			 return;
		    }
		    if (n != 0 && n != mxGetN(prhs[0])) {
			 mexErrMsgTxt("Dimension error (arg 1 and later).");
			 return;
		    }
		    n  = mxGetN(prhs[0]);
		    Q = mxGetPr(prhs[0]);
	    }
	    break;
    }
    if (nlhs > 2 || nlhs < 1) {
	    mexErrMsgTxt("Usage: [x,how] "
		         "= sqvp(H,c,l,u,x0,eq,vm)");
	    return;
    }

    switch (nlhs) {
    case 2:
    case 1:
	    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	    x = mxGetPr(plhs[0]);
	    break;
    }

    SimpleQuadraticProgram sqp(n);
    sqp.a=Q;
    sqp.b=c;
    sqp.cmin=l;
    sqp.cmax=u;
    sqp.x=x0;
    sqp.sumflag=neq;
    sqp.valmax=valmax;

    int res=sqp.run();
    int i;
    for (i=0; i<n; i++) 
      x[i]=sqp.x[i];

    if (nlhs == 2) {
	    if (res==-1) plhs[1] = mxCreateString(sqp.error);
            else plhs[1] = mxCreateString("Optimal solution found");
    }
}

