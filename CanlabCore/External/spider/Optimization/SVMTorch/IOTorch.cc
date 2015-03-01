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


#include "IOTorch.h"

void IOTorch::error()
{
  cout << "$ File error. Check your data.\n\n";
  exit(0);
}

void IOTorch::loadData(string file, real ***data_, real **y_,
                       int &l, int &w, bool bin, int max_load, bool target)
{  
   real **data;
   real *y;

  ifstream f;
  
  if(bin)
    f.open(file.c_str(), ios::in | ios::binary);
  else
    f.open(file.c_str());
    

  if(!f)
    error();

  if(bin)
  {
    f.read((char *)&l, sizeof(int));
    f.read((char *)&w, sizeof(int));
  }
  else
  {
    f >> l;
    f >> w;
  }

  if( (max_load > 0) && (max_load < l) )
  {
    l = max_load;
    cout << "# Loading only " << l << " examples\n";
  }
  if(target)
    w -= 1;
  
  data = new real*[l];
  for(int i = 0; i < l; i++)
    data[i] = new real[w];
  
  if(target)
    y = new real[l];
  else
    y = NULL;
  
  if(bin)
  {
    for(int i = 0; i < l; i++)
    {
      f.read((char *)data[i], sizeof(real)*w);
      if(target)
        f.read((char *)&y[i], sizeof(real));
    }
  }
  else
  {
    for(int i = 0; i < l; i++)
    {
      for(int j = 0; j < w; j++)
        f >> data[i][j];
      
      if(target)
        f >> y[i];
    }
  }

   *data_ = data;
   *y_ = y;
}

void IOTorch::loadData(string file, sreal ***data_, real **y_,
                       int &l, int &w, bool bin, int max_load, bool target)
{
  sreal **data;
  real *y;

  ifstream f;
  
  if(bin)
    f.open(file.c_str(), ios::in | ios::binary);
  else
    f.open(file.c_str());

  if(!f)
    error();

  if(bin)
  {
    f.read((char *)&l, sizeof(int));
    f.read((char *)&w, sizeof(int));
  }
  else
  {
    f >> l;
    f >> w;
  }

  if( (max_load > 0) && (max_load < l) )
  {
    l = max_load;
    cout << "# Loading only " << l << " examples.\n";
  }
  if(target)
    w -= 1;
  
  data = new sreal*[l];
  
  if(target)
    y = new real[l];
  else
    y = NULL;

  int w_on_line;
  if(bin)
  {
    for(int i = 0; i < l; i++)
    {
      f.read((char *)&w_on_line, sizeof(int));
      data[i] = new sreal[w_on_line+1];
      data[i][w_on_line].indice = -1;
      f.read((char *)data[i], sizeof(sreal)*w_on_line);
      if(target)
        f.read((char *)&y[i], sizeof(real));
    }
  }
  else
  {
    for(int i = 0; i < l; i++)
    {
      f >> w_on_line;
      data[i] = new sreal[w_on_line+1];
      data[i][w_on_line].indice = -1;
      for(int j = 0; j < w_on_line; j++)
      {
        f >> data[i][j].indice;
        f >> data[i][j].valeur;
      }

      if(target)
        f >> y[i];
    }
  }

  *data_ = data;
  *y_ = y;
}

//////////////// Save

void IOTorch::saveData(string file, real **data, real *y,
                       int l, int w, bool bin, int max_save, bool target)
{
  ofstream f;

  if(bin)
    f.open(file.c_str(), ios::out | ios::trunc | ios::binary);
  else
    f.open(file.c_str());

  if(!f)
    error();

  if(target)
    w++;
  if(bin)
  {
    f.write((char *)&l, sizeof(int));
    f.write((char *)&w, sizeof(int));
  }
  else
  {
    f << l << " ";
    f << w << endl;
  }
  if(target)
    w--;

  if( (max_save > 0) && (max_save < l) )
  {
    l = max_save;
    cout << "# Saving only " << l << " examples.\n";
  }
  
  if(bin)
  {
    for(int i = 0; i < l; i++)
    {
      f.write((char *)data[i], sizeof(real)*w);
      if(target)
        f.write((char *)&y[i], sizeof(real));
    }
  }
  else
  {
    for(int i = 0; i < l; i++)
    {
      for(int j = 0; j < w; j++)
      {
        f << data[i][j];
        f << " ";
      }
      
      if(target)
        f << y[i];
      f << endl;
    }
  }
}

void IOTorch::saveData(string file, sreal **data, real *y,
                       int l, int w, bool bin, int max_save, bool target)
{
  ofstream f;

  if(bin)
    f.open(file.c_str(), ios::out | ios::trunc | ios::binary);
  else
    f.open(file.c_str());

  if(!f)
    error();

  if(target)
    w++;
  if(bin)
  {
    f.write((char *)&l, sizeof(int));
    f.write((char *)&w, sizeof(int));
  }
  else
  {
    f << l << " ";
    f << w << endl;
  }
  if(target)
    w--;

  if( (max_save > 0) && (max_save < l) )
  {
    l = max_save;
    cout << "# Saving only " << l << " examples.\n";
  }

  int w_on_line;
  if(bin)
  {
    for(int i = 0; i < l; i++)
    {
      w_on_line = sparseLineLength(data[i]);
      f.write((char *)&w_on_line, sizeof(int));
      f.write((char *)data[i], sizeof(sreal)*w_on_line);
      if(target)
        f.write((char *)&y[i], sizeof(real));
    }
  }
  else
  {
    for(int i = 0; i < l; i++)
    {
      w_on_line = sparseLineLength(data[i]);
      f << w_on_line << " ";
      for(int j = 0; j < w_on_line; j++)
      {
        f << data[i][j].indice << " ";
        f << data[i][j].valeur << " ";
      }

      if(target)
        f << y[i];
      f << endl;
    }
  }
}

