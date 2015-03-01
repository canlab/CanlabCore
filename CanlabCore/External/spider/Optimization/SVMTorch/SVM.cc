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


#include "SVM.h"

SVM::SVM()
{
  sparse_mode = false;
  unshrink_mode = false;
  n_iter_min_to_shrink = 100;

#ifdef USEDOUBLE
  eps_shrink = 1E-9;
  eps_bornes = 1E-12;
#else
  eps_shrink = 1E-4;
  eps_bornes = 1E-4;
#endif

  eps_fin = 0.01;

  deja_shrink = false;

  n_max_unshrink = 1;
}

SVM::~SVM()
{
}

void SVM::setKernel(Kernel *kernel_)
{
  kernel = kernel_;
}

void SVM::setOption(const string& optionname, const string& value)
{
  if(optionname == "sparse_mode")
  {
    if(value == "1")
      sparse_mode = true;
    else
      sparse_mode = false;
  }
  else if(optionname == "unshrink_mode")
  {
    if(value == "1")
      unshrink_mode = true;
    else
      unshrink_mode = false;
  }
  else if(optionname == "n_iter_min_to_shrink")
    n_iter_min_to_shrink = atoi(value.c_str());
  else if(optionname == "n_max_unshrink")
    n_max_unshrink = atoi(value.c_str());
  else if(optionname == "eps_shrink")
    eps_shrink = atof(value.c_str());
  else if(optionname == "eps_fin")
    eps_fin = atof(value.c_str());
  else if(optionname == "eps_bornes")
    eps_bornes = atof(value.c_str());
  else
  {
    cout << "\n$ Error : unknow option >> " << optionname << " << \n\n";
    exit(0);
  }
}

bool SVM::selectVariables(int &i, int &j)
{
  real gmax_i = -INF;
  real gmin_j =  INF;
  int i_ = -1;
  int j_ = -1;

  for(int it = 0; it < n_active_var; it++)   
  {
    int t = active_var[it];

    if(y[t] > 0)
    {
      if(isNotDown(t))
      {
        if(grad[t] > gmax_i)
        {
          gmax_i = grad[t];
          i_ = t;
        }
      }

      if(isNotUp(t))
      {
        if(grad[t] < gmin_j)
        {
          gmin_j = grad[t];
          j_ = t;
        }
      }            
    }
    else
    {
      if(isNotUp(t))
      {
        if(-grad[t] > gmax_i)
        {
          gmax_i = -grad[t];
          i_ = t;
        }
      }

      if(isNotDown(t))
      {
        if(-grad[t] < gmin_j)
        {
          gmin_j = -grad[t];
          j_ = t;
        }
      }            
    }
  }

  current_error =  gmax_i - gmin_j;

  if(current_error < eps_fin)
    return(true);
  
  if( (i_ == -1) || (j_ == -1) )
    return(true);

  i = i_;
  j = j_;

  return(false);
}

bool SVM::bCompute()
{
  real sum = 0;
  int n_ = 0;
  for(int it = 0; it < n_active_var; it++)
  {
    int t = active_var[it];
    if( isNotUp(t) && isNotDown(t) )
    {
      sum += y[t]*grad[t];
      n_++;
    }
  }
  
  if(n_)
  {
    b = -sum/(real)n_;
    return(true);
  }
  else
    return(false);
}

// Renvoie le nb de var susceptibles d'etre shrinkee
int SVM::checkShrinking(real bmin, real bmax)
{
  real bb = (bmin+bmax)/2.;

  n_active_var_new = 0;
  for(int it = 0; it < n_active_var; it++)
  {
    int t = active_var[it];
    bool garde = true;

    if(isNotDown(t) && isNotUp(t))
      not_at_bound_at_iter[t] = iter;
    else
    {
      if(isNotUp(t)) // Donc elle est en bas.
      {
        if(grad[t] + y[t]*bb < eps_shrink)
          not_at_bound_at_iter[t] = iter;
        else
        {
          if( (iter - not_at_bound_at_iter[t]) > n_iter_min_to_shrink)
            garde = false;
        }
      }
      else
      {
        if(grad[t] + y[t]*bb > -eps_shrink)
          not_at_bound_at_iter[t] = iter;
        else
        {
          if( (iter - not_at_bound_at_iter[t]) > n_iter_min_to_shrink)
            garde = false;
        }
      }      
    }

    if(garde)
      active_var_new[n_active_var_new++] = t;
  }

  return(n_active_var-n_active_var_new);
}

void SVM::shrink()
{
  n_active_var = n_active_var_new;
  int *ptr_sav = active_var;
  active_var = active_var_new;
  active_var_new = ptr_sav;
  deja_shrink = true;
}

void SVM::unShrink()
{
  for(int i = 0; i < l; i++)
    active_var[i] = i;

  n_active_var = l;
  deja_shrink = false;

  if(++n_unshrink == n_max_unshrink)
  {
    unshrink_mode = false;
    n_iter_min_to_shrink = 666666666;
    cout << "shrinking and unshrinking desactived...";
  }
}

void SVM::solve()
{
  int xi, xj;
  int n_to_shrink = 0;

#ifdef I_WANT_TIME
  long t_start = getRuntime();
#endif

  n_unshrink = 0;
  b = 0;

  iter = 0;
  while(1)
  {
    if(selectVariables(xi, xj))
    {
      if(unshrink_mode)
      {
        cout << "# Unshrink...";
        cout.flush();
        unShrink();
        if(selectVariables(xi, xj))
        {
          cout << "finished.\n";
          break;
        }
        else
          cout << "restart.\n";
      }
      else
        break;
    }

//    bCompute();
//    cout << -y[xi]*grad[xi] << " < " << b << " < " << -y[xj]*grad[xj] << endl;
    if(iter >= n_iter_min_to_shrink)
      n_to_shrink = checkShrinking(-y[xi]*grad[xi], -y[xj]*grad[xj]);

    k_xi = cache->adresseCache(xi);
    k_xj = cache->adresseCache(xj);

    old_alpha_xi = alpha[xi];
    old_alpha_xj = alpha[xj];

    analyticSolve(xi, xj);

    real delta_alpha_xi = alpha[xi] - old_alpha_xi;
    real delta_alpha_xj = alpha[xj] - old_alpha_xj;

    if(deja_shrink && !unshrink_mode)
    {
      for(int t = 0; t < n_active_var; t++)
      {
        int it = active_var[t];
        grad[it] += k_xi[it]*delta_alpha_xi + k_xj[it]*delta_alpha_xj;
      }
    }
    else
    {
      for(int t = 0; t < l; t++)
        grad[t] += k_xi[t]*delta_alpha_xi + k_xj[t]*delta_alpha_xj;
    }

    iter++;
    if(! (iter % 1000) )
    {
      // Pour ne pas effrayer le neophite.
      if(current_error < 0)
        current_error = 0;
      cout << "  + Iteration " << setw(10) << iter << "\n";
      cout << "   --> Current error    = " << setw(10) << current_error << "\n";
      cout << "   --> Active variables = " << setw(10) << n_active_var  << "\n";
      cout.flush();
    }

    /////////////// Shhhhhrinnnk

    if(!(iter % n_iter_min_to_shrink))
    {
      if( (n_to_shrink > n_active_var/10) && (n_active_var-n_to_shrink > 100) )
        shrink();
    }
  }

  // Pour ne pas effrayer le neophite.
  if(current_error < 0)
    current_error = 0;
  cout << "  + Iteration " << setw(10) << iter << "\n";
  cout << "   --> Current error    = " << setw(10) << current_error << "\n";
  cout << "   --> Active variables = " << setw(10) << n_active_var  << "\n";
  cout.flush();

#ifdef I_WANT_TIME
  long t_end = getRuntime();
  cout << "# Time in CPU-seconds = " << (double)(t_end-t_start)/((double)CLK_TCK) << endl;
#endif

}

void SVM::updateStatus(int i)
{
  if(alpha[i] < Cx[i] - eps_bornes)
    status_alpha[i] = 1;
  else
    status_alpha[i] = 0;

  if(alpha[i] > eps_bornes)
    status_alpha[i] |= 2;
}
