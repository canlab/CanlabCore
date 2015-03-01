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


#ifndef STANDARDSVM_INC
#define STANDARDSVM_INC

#include "SVM.h"
#include "UserKernel.h"

class StandardSVM : public SVM
{
  private:

  real C;
  void prepareToLauch();
  void atomiseAll();
  bool svm_is_ready;
  bool support_vectors_is_ready;
  int lm;
  real eps_regression;
  real cache_size_meg;

  public:

  int n_support_vectors;
  int n_support_vectors_bound;
  int *support_vectors;
  real *sv_alpha;

  bool regression_mode;

  StandardSVM();

  virtual void analyticSolve(int xi, int xj);
  virtual void setOption(const string& optionname, const string& value);

  void train(real **data_, sreal **sdata_, real *target_, int l_, int n_input_dim_);
  real use(real *x);
  real use(sreal *x);

  void save(string file, string comment);
  void load(string file);

  ~StandardSVM();

};

#endif
