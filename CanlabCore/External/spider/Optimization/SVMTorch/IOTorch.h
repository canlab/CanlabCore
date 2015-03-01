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


#ifndef IOTORCH_INC
#define IOTORCH_INC

#include "general.h"

class IOTorch
{
  private:
  void error();

  public:
  void loadData(string file,  real ***data_, real **y_, int &l, int &w, bool bin, int max_load, bool target=true);
  void loadData(string file, sreal ***data_, real **y_, int &l, int &w, bool bin, int max_load, bool target=true);
  void saveData(string file,  real  **data,  real  *y,  int  l, int  w, bool bin, int max_save, bool target=true);
  void saveData(string file, sreal  **data,  real  *y,  int  l, int  w, bool bin, int max_save, bool target=true);
};

#endif
