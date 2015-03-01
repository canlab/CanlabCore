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


#include "general.h"
#include "StandardSVM.h"
#include "Kernel.h"
#include "UserKernel.h"
#include "IOTorch.h"

typedef struct
{
    bool regression_mode;
    bool multi_mode;
    bool sparse_mode;
    bool bin_mode;
    int  load;
    bool norm;
    bool out_ascii;
    bool out_binary;
    string out_file_a;
    string out_file_b;
    bool desired;
    int n_class;
} parametres;

void help(char **argv)
{
  cout << endl;
  cout << "==============================================================================================" << endl;
  cout << "* The Trebolloc Company presents:" << endl;
  cout << "* SVMTest" << endl;
  cout << "* a part of " << VERSION << endl;
  cout << "* \n* usage: " << argv[0] << " [options] <model file> <test file>" << endl;
  cout << "* \n* \n";
  cout << "* Arguments:" << endl;
  cout << "*   <model file> -> file to store the model" << endl;
  cout << "*   <test file>  -> file with test data" << endl;
  cout << "* General options:" << endl;
  cout << "*   -help        -> this little help" << endl;
  cout << "*   -multi       -> multiclass mode" << endl;
  cout << "*   -norm        -> normalize in multiclass mode" << endl;
  cout << "* Input option:" << endl;
  cout << "*   -sparse      -> sparse data" << endl;
  cout << "*   -bin         -> binary data" << endl;
  cout << "*   -load <int>  -> load only <int> test examples" << endl;
  cout << "*   -no          -> no desired-output in test_file" << endl;
  cout << "* Output options:" << endl;
  cout << "*   -oa <file>   -> write in ASCII the SVM output into <file>" << endl;
  cout << "*   -ob <file>   -> write in binary the SVM output into <file>" << endl;
  cout << "*" << endl;
  cout << "===========>> SVMTorch is (c) Ronan Collobert 2001 (IDIAP & DIRO U. de Montreal) <<===========" << endl;
  cout << "\n";
  exit(0);
}

void scan_cmd(int argc, char **argv, parametres *params)
{  
  int i = 1;
  int erreur;
  int maxa;
  
  if(argc < 3)
    help(argv);
  
  maxa = argc-2;
  while(i < maxa)
  {
    erreur = 1;
    
    if(!strcmp(argv[i], "-help"))
      help(argv);

    if(!strcmp(argv[i], "-multi"))
    {
      params->multi_mode = true;
      erreur = 0;
    }    

    if(!strcmp(argv[i], "-sparse"))
    {
      params->sparse_mode = true;
      erreur = 0;
    }

    if(!strcmp(argv[i], "-norm"))
    {
      params->norm = true;
      erreur = 0;
    }

    if(!strcmp(argv[i], "-oa"))
    {
      i++;
      if(i < maxa)
        params->out_file_a = argv[i];
      else
        help(argv);

      params->out_ascii = true;      
      erreur = 0;
    }

    if(!strcmp(argv[i], "-ob"))
    {
      i++;
      if(i < maxa)
        params->out_file_b = argv[i];
      else
        help(argv);

      params->out_binary = true;      
      erreur = 0;
    }

    if(!strcmp(argv[i], "-bin"))
    {
      params->bin_mode = true;
      erreur = 0;
    }

    if(!strcmp(argv[i], "-no"))
    {
      params->desired = false;
      erreur = 0;
    }

    if(!strcmp(argv[i], "-load"))
    {
      i++;
      if(i < maxa)
        params->load = atoi(argv[i]);
      else
        help(argv);
      
      erreur = 0;
    }

    i++;    
    if(erreur)
      help(argv);
  }
}

real rabouleIllico(string file_svm, real **data, sreal **sdata, real *y_pred, int l, parametres *params)
{
  StandardSVM emilie;
  
  emilie.load(file_svm);

  if(emilie.regression_mode)
  {
    params->regression_mode = true;
    cout << "# Regression mode" << endl;
  }
  else
    cout << "# Classification mode" << endl;


  int n_aff_prochain = 0;
  int n_deja_aff = 0;
  cout << "                    ";
  cout << "[________Test________]" << endl;
  cout << "                    [";
  cout.flush();
 
  if(params->sparse_mode)
  {
    for(int i = 0; i < l; i++)
    {
      y_pred[i] = emilie.use(sdata[i]);

      if(i >= n_aff_prochain)
      {
        n_aff_prochain = ++n_deja_aff * l / 20;
        cout << "#";
        cout.flush();
      }
    }
  }
  else
  {
    for(int i = 0; i < l; i++)
    {
      y_pred[i] = emilie.use(data[i]);

      if(i >= n_aff_prochain)
      {
        n_aff_prochain = ++n_deja_aff * l / 20;
        cout << "#";
        cout.flush();
      }
    }
  }
  cout << "]" << endl;

  if(params->multi_mode)
  {
    real sum = 0;
    for(int i = 0; i < emilie.n_support_vectors; i++)
      sum += fabs(emilie.sv_alpha[i]);

    return(sum);
  }

  return(0.);
}


int main(int argc, char **argv)
{
  int *multi_missclassified = NULL;
  int *multi_missclassified_pos = NULL;

  parametres params;

  params.multi_mode = false;
  params.sparse_mode = false;
  params.bin_mode = false;    
  params.load = -1;
  params.norm = false;
  params.out_ascii = false;
  params.out_binary = false;
  params.out_file_a = "zarma";
  params.out_file_b = "ether";
  params.desired = true;
  params.regression_mode = false;
  params.n_class = -1;
  scan_cmd(argc, argv, &params);

  if(!params.desired)
  {
    if(!params.out_ascii && !params.out_binary)
    {
      cout << "$ What do you want boy ??" << endl;
      exit(0);
    }
  }

  // Let's go boy //////////////////////////////////
  real **data = NULL;
  sreal **sdata = NULL;
  real *y = NULL;
  int l, c;

  cout << "# Loading data :" << endl;
  IOTorch melanie;
  if(params.sparse_mode)
    melanie.loadData(argv[argc-1], &sdata, &y, l, c, params.bin_mode, params.load, params.desired);
  else
    melanie.loadData(argv[argc-1], &data, &y, l, c, params.bin_mode, params.load, params.desired);

  cout << "  " << l << " test examples" << endl;
  cout << "  input dimension is " << c << endl;

  cout << "# Starting test..." << endl;
  real *y_pred = new real[l];
  if(params.multi_mode)
  {
    // Postulat J4 : "A tout probleme il y a une solution boeuf" ////
    cout << "# Scanning for classes...";
    cout.flush();
    int n_class = 0;

    while(1)
    {
      string file = argv[argc-2];
      {
        char julie[10];
        sprintf(julie, ".%d", n_class);
        string sjulie = julie;
        file += sjulie;
      }
      
      ifstream f(file.c_str());
      if(!f)
        break;
      else
        n_class++;
    }
    params.n_class = n_class;
    cout << n_class << " found" << endl;
    /////////////

    int *class_pred = new int[l];
    for(int i = 0; i < l; i++)
    {
      y_pred[i] = -INF;
      class_pred[i] = -1;
    }

    //////////////////////////////////////////////
    if(params.desired)
    {
      multi_missclassified = new int[n_class];
      multi_missclassified_pos = new int[n_class];
    }
    //////////////////////////////////////////////

    real *y_temp = new real[l];

    for(int cl = 0; cl < n_class; cl++)
    {
      cout << "\n# Testing class " << cl << endl;
      string file = argv[argc-2];
      {
        char julie[10];
        sprintf(julie, ".%d", cl);
        string sjulie = julie;
        file += sjulie;
      }
      
      real norm = rabouleIllico(file, data, sdata, y_temp, l, &params);
      if(params.norm)
      {
        for(int i = 0; i < l; i++)
          y_temp[i] /= norm;

        cout << "# Norm = " << norm << endl;
      }

      // Bon, faut bien les compter la... /////
      if(params.desired)
      {
        multi_missclassified[cl] = 0;
        multi_missclassified_pos[cl] = 0;
        for(int i = 0; i < l; i++)
        {
          real z = ( ((int)y[i] == cl) ? 1 : -1);
          if(z*y_temp[i] <= 0)
          {
            multi_missclassified[cl]++;
            if(z < 0)
              multi_missclassified_pos[cl]++;
          }       
        }
      }
      /////////////////////////////////////////

      for(int i = 0; i < l; i++)
      {
        if(y_temp[i] > y_pred[i])
        {
          y_pred[i] = y_temp[i];
          class_pred[i] = cl;
        }
      }
    }
    for(int i = 0; i < l; i++)
      y_pred[i] = (real)class_pred[i];
    delete[] y_temp;
    delete[] class_pred;
  }
  else
    rabouleIllico(argv[argc-2], data, sdata, y_pred, l, &params);


  /// Bon. Les calculs divers et debiles.

  // Pour la regression
  if(params.regression_mode && params.desired && params.desired)
  {
    real maxe = 0;
    real mae = 0;
    real mse = 0;
    real msse = 0;
    
    for(int i = 0; i < l; i++)
    {
      real z = y[i] - y_pred[i];
      z = ( (z < 0) ? -z : z);
      if(z < 0)
        z = -z;
            
      if(z > maxe)
        maxe = z;
      
      mae += z;
      z *= z;
      mse += z;
      z *= z;
      msse += z;
    }

    {
      real z = (real)l;
      mae /= z;
      mse /= z;
      msse /= z;
    }
    
    cout << endl;
    cout << "# Mean Absolute Error    : " << mae << endl;
		cout << "   -> Standard Deviation : " << sqrt(mse-mae*mae) << endl;
    cout << endl;
    cout << "# Mean Squared Error     : " << mse << endl;
    cout << "   -> Standard Deviation : " << sqrt(msse-mse*mse) << endl;
    cout << endl;
    cout << "# Maximal Error          : " << maxe << endl;
  }

  // Pour la classification a deux classes.
  if(!params.regression_mode && !params.multi_mode && params.desired)
  {
    int missclassified = 0;
    int missclassified_pos = 0;
    for(int i = 0; i < l; i++)
    {
      if(y[i]*y_pred[i] <= 0)
      {
        missclassified++;
        if(y[i] < 0)
          missclassified_pos++;
      }
    }
    
    real z = 100.0 * ((real)missclassified) / ((real)l);
    cout << endl;
		cout << "# Number of missclassified          : " << missclassified << " [" << setprecision(4) << z << "%]" << endl;
		cout << "     -> False positives             : " << missclassified_pos << endl;
		cout << "     -> False negatives             : " << missclassified-missclassified_pos << endl;
    z = 100.0 * ((real)(l-missclassified)) / ((real)l);
    cout << "# Number of correct classifications : " << l-missclassified << " ["  << setprecision(4) <<  z << "%]" << endl;
  }

  // Pour la classification multiclasse.
  if(!params.regression_mode && params.multi_mode && params.desired)
  {
    int missclassified = 0;
    for(int i = 0; i < l; i++)
    {
      if(y[i] != y_pred[i])
        missclassified++;
    }

    cout << endl;
    cout << "# Number of missclassified for each class:" << endl;
    for(int i = 0; i < params.n_class; i++)
    {
      real z = 100.0 * ((real)multi_missclassified[i]) / ((real)l);
      cout << "  + For class " << setw(4) << i << "                  : " << multi_missclassified[i] << " [" << setprecision(4) << z << "%]" << endl;
      cout << "     -> False positives             : " << multi_missclassified_pos[i] << endl;
      cout << "     -> False negatives             : " << multi_missclassified[i]-multi_missclassified_pos[i] << endl;
      cout << endl;
    }

    cout << endl;
    cout << "# With multiclass: " << missclassified << " missclassified [" <<  setprecision(4) << (100. * (real)missclassified / (real)l);
    cout << "%] on " << l << " examples" << endl;    
    delete[] multi_missclassified;
    delete[] multi_missclassified_pos;
  }  
  cout << endl;

  // Sauvegardes diverses...

  // Ascii
  if(params.out_ascii)
  {
    cout << "# Writing results in ASCII" << endl;
    ofstream f( (params.out_file_a).c_str() );
    for(int i = 0; i < l; i++)
      f << y_pred[i] << endl;
  }

  if(params.out_binary)
  {
    cout << "# Writing results in binary" << endl;
    ofstream f( (params.out_file_b).c_str() );
    f.write((char *)y_pred, sizeof(real)*l);
  }

  cout << endl;

  return(0);
}


