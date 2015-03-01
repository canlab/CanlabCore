#ifndef VECTORUTIL
#define VECTORUTIL
#include <stdlib.h>

/**
 * @author froehlic
 *
 * some useful vector mathematics for simple array representations of vectos
 
	/**
   * multiply each element of @param v by @param s
   */
  double* svmult (double s, double* v, int nv) {
    double* res = new double[nv];
    for (int i = 0; i < nv; i++) {
        res[i] = v[i]*s;
    }
    return  res;
  }
  /**
   * divide each element of @param v by @param s
   */
  double* svdiv (double s, double* v, int nv) {
    double* res = new double[nv];
    for (int i = 0; i < nv; i++) {
        res[i] = v[i]/s;
    }
    return res;
  }
  /**
   * add @param s to each element of @param v
   */
  double* svadd (double s, double* v, int nv) {
    double* res = new double[nv];
    for (int i = 0; i < nv; i++) {
        res[i] = v[i] + s;
    }
    return res;
  }
  /**
   * @return @param v1 + @parm v2
   */
  double* vvadd (double* v1, double* v2, int nv) {
    double* res = new double[nv];
    for (int i = 0; i < nv; i++) {
        res[i] = v1[i] + v2[i];
    }
    return res;
  }
  /**
   * @return @param v1 - @param v2 
   */
  double* vvsub (double* v1, double* v2, int nv) {
    double* res = new double[nv];
    for (int i = 0; i < nv; i++) {
        res[i] = v1[i] - v2[i];
    }
    return res;
  }
  /**
   * @return scalar product of @param v1 and @param v2
   */
  double scalar(double* v1, double* v2, int nv) {
    double res = 0;
    for (int i = 0; i < nv; i++) {
        res += v1[i] * v2[i];
    }
    return res;
  }
/**
 * clone @param values
 *  */
 double* clone(double* values, int nv){
  	double* result = new double[nv];
  	for(int i = 0; i < nv; i++)
  		result[i] = values[i];
  	return result;
  }
  
#endif
