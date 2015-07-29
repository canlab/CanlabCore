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


#include "Kernel.h"

Kernel::Kernel(SVM &svm_) : svm(svm_)
{
  svm.setKernel(this);
}

real Kernel::dot(int i, int j)
{
  if(svm.sparse_mode)
  {
    real sum = 0;
    sreal *ptr_i = svm.sdata[i];
    sreal *ptr_j = svm.sdata[j];
    while( (ptr_i->indice != -1) && (ptr_j->indice != -1) )
    {
      if(ptr_i->indice == ptr_j->indice)
      {
        sum += ptr_i->valeur * ptr_j->valeur;
        ptr_i++;
        ptr_j++;
      }
      else
      {
        if(ptr_i->indice > ptr_j->indice)
          ptr_j++;
        else
          ptr_i++;
      }
    }
      
    return(sum);      
  }
  else
  {
    real sum = 0;
    real *ptr_i = svm.data[i];
    real *ptr_j = svm.data[j];
    for(int t = 0; t < svm.n_input_dim; t++)
      sum += ptr_i[t] * ptr_j[t];

    return(sum);
  }
}

void Kernel::init()
{
}

real Kernel::dotStandard(real *x, int j)
{
  real sum = 0;
  real *ptr_i = x;
  real *ptr_j = svm.data[j];
  for(int t = 0; t < svm.n_input_dim; t++)
    sum += ptr_i[t] * ptr_j[t];

  return(sum);
}

real Kernel::dotSparse(sreal *x, int j)
{
  real sum = 0;
  sreal *ptr_i = x;
  sreal *ptr_j = svm.sdata[j];
  while( (ptr_i->indice != -1) && (ptr_j->indice != -1) )
  {
    if(ptr_i->indice == ptr_j->indice)
    {
      sum += ptr_i->valeur * ptr_j->valeur;
      ptr_i++;
      ptr_j++;
    }
    else
    {
      if(ptr_i->indice > ptr_j->indice)
        ptr_j++;
      else
        ptr_i++;
    }
  }

  return(sum);
}

Kernel::~Kernel()
{

}

////
DotKernel::DotKernel(SVM &svm_) : Kernel(svm_)
{
  id = 0;
}

real DotKernel::evalue(int i, int j)
{
  return(dot(i, j));
}

real DotKernel::evalue(real *x, int j)
{
  return(dotStandard(x, j));
}

real DotKernel::evalue(sreal *x, int j)
{
  return(dotSparse(x, j));
}

////
PolynomialKernel::PolynomialKernel(SVM &svm_) : Kernel(svm_)
{
  id = 1;
  setParametres();
}

void PolynomialKernel::setParametres(int d_, real s_, real r_)
{
  d = d_;
  s = s_;
  r = r_;
}

void PolynomialKernel::getParametres(int &d_, real &s_, real &r_)
{
  d_ = d;
  s_ = s;
  r_ = r;
}

real PolynomialKernel::evalue(int i, int j)
{
  real z = s*dot(i, j) + r;

  // la fonction pow rame a donf
  real julie = z;
  for(int t = 1; t < d; t++)
    julie *= z;
    
  return(julie);
}

real PolynomialKernel::evalue(real *x, int j)
{
  real z = s*dotStandard(x, j) + r;

  // la fonction pow rame a donf
  real julie = z;
  for(int t = 1; t < d; t++)
    julie *= z;
    
  return(julie);
}

real PolynomialKernel::evalue(sreal *x, int j)
{
  real z = s*dotSparse(x, j) + r;

  // la fonction pow rame a donf
  real julie = z;
  for(int t = 1; t < d; t++)
    julie *= z;
    
  return(julie);
}

////
GaussianKernel::GaussianKernel(SVM &svm_) : Kernel(svm_)
{
  id = 2;
  setParametres();
  precalc_alloc = false;
}

void GaussianKernel::init()
{
  precalc = new real[svm.n_train_examples];
  cout << "# Precalculating...";
  cout.flush();
  for(int i = 0; i < svm.n_train_examples; i++)
    precalc[i] = dot(i, i);
  cout << "OK\n";

  precalc_alloc = true;
}

void GaussianKernel::setParametres(real std)
{
//  g = 1./(2.*std*std);
  g = 1./(std*std);
}

void GaussianKernel::getParametres(real &std)
{
//  std = 1./(sqrt(2*g));
  std = 1./(sqrt(g));
}

real GaussianKernel::evalue(int i, int j)
{    
  return(exp(g*( 2.*dot(i,j)-precalc[i]-precalc[j] )));
}

real GaussianKernel::evalue(real *x, int j)
{
  real sum = 0;
  real *y = svm.data[j];

  for(int t = 0; t < svm.n_input_dim; t++)
  {
    real z = x[t] - y[t];
    sum -= z*z;
  }

  return(exp(g*sum));
}

real GaussianKernel::evalue(sreal *x, int j)
{

  real sum = 0;
  sreal *ptr_i = x;
  sreal *ptr_j = svm.sdata[j];
  while( (ptr_i->indice != -1) && (ptr_j->indice != -1) )
  {
    if(ptr_i->indice == ptr_j->indice)
    {
      real z = ptr_i->valeur - ptr_j->valeur;
      sum += z*z;
      ptr_i++;
      ptr_j++;
    }
    else
    {
      if(ptr_i->indice > ptr_j->indice)
      {
        sum += ptr_j->valeur * ptr_j ->valeur;
        ptr_j++;
      }
      else
      {
        sum += ptr_i->valeur * ptr_i ->valeur;
        ptr_i++;
      }
    }
  }

  while(ptr_i->indice != -1)
  {
    sum += ptr_i->valeur * ptr_i ->valeur;
    ptr_i++;
  }

  while(ptr_j->indice != -1)
  {
    sum += ptr_j->valeur * ptr_j ->valeur;
    ptr_j++;
  }

  return(exp(-g*sum));
}

GaussianKernel::~GaussianKernel()
{
  if(precalc_alloc)
    delete[] precalc;
}


////
SigmoidKernel::SigmoidKernel(SVM &svm_) : Kernel(svm_)
{
  id = 3;
  setParametres();
}

void SigmoidKernel::setParametres(real s_, real r_)
{
  s = s_;
  r = r_;
}

void SigmoidKernel::getParametres(real &s_, real &r_)
{
  s_ = s;
  r_ = r;
}

real SigmoidKernel::evalue(int i, int j)
{
  return(tanh(s*dot(i, j) + r));
}

real SigmoidKernel::evalue(real *x, int j)
{
  return(tanh(s*dotStandard(x, j) + r));
}

real SigmoidKernel::evalue(sreal *x, int j)
{
  return(tanh(s*dotSparse(x, j) + r));
}
