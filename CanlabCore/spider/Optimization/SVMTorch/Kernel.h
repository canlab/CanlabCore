//===////////////////////////////////////////////////////////////////////
//                                                                     //
//   SVMTorch II V1.77 [Ennastie]                                      //
//                                                                     //
//   Author: Ronan Collobert                                           //
//   Date: 20.11.2001                                                  //
//                                                                     //
//   Copyright (c) 2001  IDIAP & DIRO U. de Montreal                   //
//                       All rights reserved                           //
//                                                                     //
//   This software is available for non-commercial use only. It must   //
//   not be modified and distributed without prior permission of the   //
//   author. The author is not responsible for implications from the   //
//   use of this software.                                             //
//                                                                     //
////////////////////////////////////////////////////////////////////===//


#ifndef INC_KERNEL
#define INC_KERNEL

#include "general.h"
#include "SVM.h"

class SVM;

class Kernel
{
  protected:
  SVM &svm;

  real dot(int i, int j);
  real dotStandard(real *x, int j);
  real dotSparse(sreal *x, int j);

  public:
  int id;
  Kernel(SVM &svm_);
  virtual void init();
  virtual real evalue(int i, int j)= 0;
  virtual real evalue(real *x, int j)= 0;
  virtual real evalue(sreal *x, int j)= 0;
  virtual ~Kernel();
};

class DotKernel : public Kernel
{
  public:
  DotKernel(SVM &svm_);
  virtual real evalue(int i, int j);
  virtual real evalue(real *x, int j);
  virtual real evalue(sreal *x, int j);
};


class PolynomialKernel : public Kernel
{
  private:
  int d;
  real s, r;

  public:
  PolynomialKernel(SVM &svm_);
  void setParametres(int d_=2, real s_=1., real r_=1.);
  void getParametres(int &d_, real &s_, real &r_);
  virtual real evalue(int i, int j);
  virtual real evalue(real *x, int j);
  virtual real evalue(sreal *x, int j);
};

class GaussianKernel : public Kernel
{
  private:
  real g;
  real *precalc;
  bool precalc_alloc;

  public:
  GaussianKernel(SVM &svm_);
  virtual void init();
  void setParametres(real std=10.);
  void getParametres(real &std);
  virtual real evalue(int i, int j);
  virtual real evalue(real *x, int j);
  virtual real evalue(sreal *x, int j);
  virtual ~GaussianKernel();
};

class SigmoidKernel : public Kernel
{
  private:
  real s, r;

  public:
  SigmoidKernel(SVM &svm_);
  void setParametres(real s_=1., real r_=1.);
  void getParametres(real &s_, real &r_);
  virtual real evalue(int i, int j);
  virtual real evalue(real *x, int j);
  virtual real evalue(sreal *x, int j);
};

#endif
