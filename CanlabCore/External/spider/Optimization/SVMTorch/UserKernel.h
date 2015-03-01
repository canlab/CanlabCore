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


#ifndef USERKERNEL_INC
#define USERKERNEL_INC

#include "Kernel.h"

class UserKernel : public Kernel
{
  private:
  real g;
  real *precalc;
  bool precalc_alloc;

  public:
  UserKernel(SVM &svm_);
  virtual void init();
  void setParametres(string u);
  void getParametres(string &u);
  virtual real evalue(int i, int j);
  virtual real evalue(real *x, int j);
  virtual real evalue(sreal *x, int j);
  virtual ~UserKernel();
};

#endif
