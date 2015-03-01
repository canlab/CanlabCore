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


#include "StandardSVM.h"

StandardSVM::StandardSVM() : SVM()
{
  C = 100;
  svm_is_ready = false;
  support_vectors_is_ready = false;
  n_support_vectors = 0;
  n_support_vectors_bound = 0;
  eps_regression = 0.7;
  cache_size_meg = 50;
  regression_mode = false;
}

StandardSVM::~StandardSVM()
{
  atomiseAll();
}

void StandardSVM::prepareToLauch()
{
  cout << "# Squatting memory\n";

  n_active_var = l;
  active_var = new int[l];
  active_var_new = new int[l];
  grad = new real[l];
  not_at_bound_at_iter = new int[l];
  alpha = new real[l];
  status_alpha = new char[l];
  Cx = new real[l];

  for(int i = 0; i < l; i++)
  {
    active_var[i] = i;
    alpha[i] = 0;
    status_alpha[i] = 1;
    not_at_bound_at_iter[i] = 0;
    Cx[i] = C;
  }

  if(regression_mode)
  {
    y = new real[l];
    for(int i = 0 ; i < lm; i++)
    {
      grad[i] =  target[i   ] + eps_regression;
      y[i] =  1;
    }
    for(int i = lm; i < l ; i++)
    {
      grad[i] = -target[i-lm] + eps_regression;
      y[i] = -1;
    }
  }
  else
  {
    for(int i = 0; i < l; i++)
      grad[i] = -1;
  }

  kernel->init();

  if(regression_mode)
    cache = new CacheRegression(*this, cache_size_meg);
  else
    cache = new CacheClassification(*this, cache_size_meg);

  svm_is_ready = true;
}

void StandardSVM::atomiseAll()
{
  if(!svm_is_ready)
    return;

  delete[] active_var;
  delete[] active_var_new;
  delete[] grad;
  delete[] not_at_bound_at_iter;
  delete[] alpha;
  delete[] status_alpha;
  delete[] Cx;
  if(regression_mode)
    delete[] y;

  delete cache;

  svm_is_ready = false;

  if(support_vectors_is_ready)
  {
    delete[] support_vectors;
    delete[] sv_alpha;
    n_support_vectors = 0;
    support_vectors_is_ready = false;
  }
}

void StandardSVM::train(real **data_, sreal **sdata_, real *target_, int l_, int n_input_dim_)
{
  atomiseAll();

  if(regression_mode)
  {
    cout << "# Regression mode\n";
    l = 2*l_;
    lm = l_;
  }
  else
  {
    cout << "# Classification mode\n";
    l = l_;
    lm = l_;
  }

  n_input_dim = n_input_dim_;
  n_train_examples = l_;

  if(sparse_mode)
    sdata = sdata_;
  else
    data = data_;

  target = target_;
  if(!regression_mode)
    y = target;


  prepareToLauch();

  cout << "# System loaded\n";
  solve();
  cout << "# System thermonuclearised\n";

  if(!bCompute())
  {
    cout << "! Warning : b is not unique. It's probably wrong.\n";
    cout << "! I think you are using silly parameters.\n";
  }

  if(regression_mode)
    b = -b;

  // La je prie pour que l'utilisateur normal utilise
  // un processeur deterministe ///
  n_support_vectors = 0;
  n_support_vectors_bound = 0;
  for(int i = 0; i < l; i++)
  {
    if(alpha[i] > eps_bornes)
    {
      if(alpha[i] > C - eps_bornes)
        n_support_vectors_bound++;

      n_support_vectors++;
    }
  }
  support_vectors = new int[n_support_vectors];
  sv_alpha = new real[n_support_vectors];

  n_support_vectors = 0;
  for(int i = 0; i < l; i++)
  {
    if(alpha[i] > eps_bornes)
    {
      support_vectors[n_support_vectors] = i;
      if(regression_mode)
      {
        if(i < n_train_examples)
          sv_alpha[n_support_vectors++] = -alpha[i];
        else
          sv_alpha[n_support_vectors++] =  alpha[i];
      }
      else
        sv_alpha[n_support_vectors++] = y[i]*alpha[i];
    }
  }

  /////

  cout << "# " << n_support_vectors << " support vectors\n";
  cout << "# With " << n_support_vectors_bound << " support vectors at C\n";
  support_vectors_is_ready = true;
}


void StandardSVM::setOption(const string& optionname, const string& value)
{
  if(optionname == "C")
    C = atof(value.c_str());
  else if(optionname == "eps_regression")
    eps_regression = atof(value.c_str());
  else if(optionname == "cache_size_meg")
    cache_size_meg = atof(value.c_str());
  else if(optionname == "regression_mode")
  {
    if(value == "1")
      regression_mode = true;
    else
      regression_mode = false;
  }
  else
    SVM::setOption(optionname, value);
}


void StandardSVM::analyticSolve(int xi, int xj)
{
  real ww, H, L;

  real s = y[xi]*y[xj];
  if(s < 0)
  {
    ww = old_alpha_xi - old_alpha_xj;
    L = ((ww   > 0.0) ? ww :  0.0);
    H = ((C+ww >   C) ? C  : C+ww);
  }
  else
  {
    ww = old_alpha_xi + old_alpha_xj;
    L = ((ww-C > 0.0) ? ww-C : 0.0);
    H = ((ww   >   C) ? C    :  ww);
  }

    
  real eta = k_xi[xi] - 2.*s*k_xi[xj] + k_xj[xj];

  if(eta > 0)
  {
    real alph = old_alpha_xi + (s*grad[xj] - grad[xi])/eta;
	
    if(alph > H)
      alph = H;
    else
    {
      if(alph < L)
        alph = L;
    }
    
    alpha[xi] = alph;
    alpha[xj] -= s*(alpha[xi]-old_alpha_xi);
  }
  else
  {
    real alph = grad[xi] - s*grad[xj];
    if(alph > 0)
    {
      alpha[xi] = L;
      alpha[xj] += s*(alpha[xi]-old_alpha_xi);
    }
    else
    {
      alpha[xi] = H;
      alpha[xj] += s*(alpha[xi]-old_alpha_xi);
    }
  }

  updateStatus(xi);
  updateStatus(xj);
}

void StandardSVM::save(string file, string comment)
{
  if(!support_vectors_is_ready)
  {
    cout << "$ Nothing to save.\n\n";
    exit(0);
  }

  ofstream f(file.c_str(), ios::out | ios::trunc | ios::binary);

  if(!f)
  {
    cout << "$ File error. Arg.\n" << endl;
    exit(0);
  }

  f << "# " << VERSION << endl;
  f << "# " << comment << endl;
  f << "# " << n_support_vectors << " support vectors inside" << endl;
  f << "# With " << n_support_vectors_bound << " support vectors at C" << endl;

  f.write((char *)&regression_mode, sizeof(bool));
  f.write((char *)&sparse_mode, sizeof(bool));

  // Kernel save
  {
    real r, s, std;
    int d;
    string u;

    f.write((char *)&kernel->id, sizeof(int));
    switch(kernel->id)
    {
      case 0:
        break;
      case 1:
        ((PolynomialKernel *)kernel)->getParametres(d, s, r);
        f.write((char *)&d, sizeof(int));
        f.write((char *)&s, sizeof(real));
        f.write((char *)&r, sizeof(real));        
        break;
      case 2:
        ((GaussianKernel *)kernel)->getParametres(std);
        f.write((char *)&std, sizeof(real));
        break;
      case 3:
        ((SigmoidKernel *)kernel)->getParametres(s, r);
        f.write((char *)&s, sizeof(real));
        f.write((char *)&r, sizeof(real));
        break;
      case 4:
        ((UserKernel *)kernel)->getParametres(u);
        f << u << endl;
        break;
      default:
        cout << "$ C'est quoi ce foutu kernel, bordel d'ostie ?!?\n\n";
        exit(0);
    }
  }

  f.write((char *)&n_support_vectors, sizeof(int));
  f.write((char *)&n_input_dim, sizeof(int));
  f.write((char *)&b, sizeof(real));

  if(sparse_mode)
  {
    int n_on_line;
    for(int it = 0; it < n_support_vectors; it++)
    {
      int t = support_vectors[it];
      n_on_line = sparseLineLength(sdata[t%lm]);
      f.write((char *)&sv_alpha[it], sizeof(real));      
      f.write((char *)&n_on_line, sizeof(int));
      f.write((char *)sdata[t%lm], sizeof(sreal)*n_on_line);
    }
  }
  else
  {
    for(int it = 0; it < n_support_vectors; it++)
    {
      int t = support_vectors[it];
      f.write((char *)&sv_alpha[it], sizeof(real));
      f.write((char *)data[t%lm], sizeof(real)*n_input_dim);
    }
  }
}

void StandardSVM::load(string file)
{
  ifstream f(file.c_str(), ios::in | ios::binary);

  if(!f)
  {
    cout << "$ File error. Arg.\n" << endl;
    exit(0);
  }

  // Jarte les commentaires
  {
    char *buffer = new char[1000];
    while(f.peek() == '#')
      f.getline(buffer, 1000);
    delete[] buffer;
  }

  f.read((char *)&regression_mode, sizeof(bool));
  f.read((char *)&sparse_mode, sizeof(bool));

  // Kernel load
  {
    real r, s, std;
    int d;
    char *u = new char[1000];
    int id;

    f.read((char *)&id, sizeof(int));
    switch(id)
    {
      case 0:
        kernel = new DotKernel(*this);
        break;
      case 1:
        kernel = new PolynomialKernel(*this);
        f.read((char *)&d, sizeof(int));
        f.read((char *)&s, sizeof(real));
        f.read((char *)&r, sizeof(real));
        ((PolynomialKernel *)kernel)->setParametres(d, s, r);
        break;
      case 2:
        kernel = new GaussianKernel(*this);
        f.read((char *)&std, sizeof(real));
        ((GaussianKernel *)kernel)->setParametres(std);
        break;
      case 3:
        kernel = new SigmoidKernel(*this);
        f.read((char *)&s, sizeof(real));
        f.read((char *)&r, sizeof(real));
        ((SigmoidKernel *)kernel)->setParametres(s, r);
        break;
      case 4:
        kernel = new UserKernel(*this);
        f.getline(u, 1000);
        ((UserKernel *)kernel)->setParametres(u);
        break;
      default:
        cout << "$ C'est quoi ce foutu kernel, bordel d'ostie ?!?\n\n";
        exit(0);
    }

    delete[] u;
  }

  f.read((char *)&n_support_vectors, sizeof(int));
  support_vectors = new int[n_support_vectors];
  sv_alpha = new real[n_support_vectors];
  f.read((char *)&n_input_dim, sizeof(int));
  f.read((char *)&b, sizeof(real));

  if(sparse_mode)
  {
    sdata = new sreal*[n_support_vectors];
    int n_on_line;
    for(int t = 0; t < n_support_vectors; t++)
    {
      support_vectors[t] = t;
      f.read((char *)&sv_alpha[t], sizeof(real));      
      f.read((char *)&n_on_line, sizeof(int));
      sdata[t] = new sreal[n_on_line+1];
      sdata[t][n_on_line].indice = -1;
      f.read((char *)sdata[t], sizeof(sreal)*n_on_line);
    }
  }
  else
  {
    data = new real*[n_support_vectors];
    for(int t = 0; t < n_support_vectors; t++)
    {
      support_vectors[t] = t;
      f.read((char *)&sv_alpha[t], sizeof(real));
      data[t] = new real[n_input_dim];
      f.read((char *)data[t], sizeof(real)*n_input_dim);
    }
  }
  l = n_support_vectors;
  lm = n_support_vectors;
  n_train_examples = l;
  kernel->init();
  cout << "# " << n_support_vectors << " support vectors in the model" << endl;
  support_vectors_is_ready = true;
}

real StandardSVM::use(real *x)
{
  if(!support_vectors_is_ready)
  {
    cout << "$ Error. No support vectors inside.\n\n";
    exit(0);
  }

  real sum = 0;
  for(int it = 0; it < n_support_vectors; it++)
  {
    int t = support_vectors[it];
    sum += sv_alpha[it]*kernel->evalue(x, t%lm);
  }

  sum += b;

  return(sum);
}


real StandardSVM::use(sreal *x)
{
  if(!support_vectors_is_ready)
  {
    cout << "$ Error. No support vectors inside.\n\n";
    exit(0);
  }

  real sum = 0;
  for(int it = 0; it < n_support_vectors; it++)
  {
    int t = support_vectors[it];
    sum += sv_alpha[it]*kernel->evalue(x, t%lm);
  }

  sum += b;

  return(sum);
}
