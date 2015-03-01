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


#ifndef INC_SVM
#define INC_SVM

#include "Kernel.h"
#include "Cache.h"

class Kernel;

class SVM
{
  friend class Cache;

  private:
  int n_unshrink;
  int n_max_unshrink;

  protected:
  real *k_xi;
  real *k_xj;
  
  real old_alpha_xi;
  real old_alpha_xj;
  real current_error;

  int *active_var_new;
  int n_active_var_new;

  public:
  int l;                  // le nb de alphas
  int n_train_examples;   // c'est parlant non ?
  Kernel *kernel;
  Cache *cache;

  bool sparse_mode;
  bool deja_shrink;
  bool unshrink_mode;

  real **data;
  sreal **sdata;
  real *y;
  real *target;
  int n_input_dim;
  
  real *alpha;
  
  real *grad;
  real eps_shrink;
  real eps_fin;
  real eps_bornes;

  real b;
  
  int n_active_var;
  int *active_var;
  int *not_at_bound_at_iter;
  int iter;
  int n_iter_min_to_shrink;

  char *status_alpha;
  real *Cx;

  SVM();
  bool bCompute();
  bool selectVariables(int &i, int &j);
  int checkShrinking(real bmin, real bmax);
  void shrink();
  void unShrink();
  void solve();
  void setKernel(Kernel *kernel_);

  virtual void analyticSolve(int xi, int xj) = 0;
  virtual void setOption(const string& optionname, const string& value);

  void updateStatus(int i);  
  inline bool isNotUp(int i)   {  return(status_alpha[i] != 2); };
  inline bool isNotDown(int i) {  return(status_alpha[i] != 1); };

  virtual ~SVM();


};

#endif
