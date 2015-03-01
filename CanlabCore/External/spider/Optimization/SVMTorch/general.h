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


#ifndef INC_GENERAL
#define INC_GENERAL

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdlib>

#include <stdio.h>
#include <limits.h>

// Old systems need that to define FLT_MAX and DBL_MAX
#ifndef DBL_MAX
#include <values.h>
#endif

using namespace std;

#define INF DBL_MAX
#define USEDOUBLE
#define I_WANT_TIME

#ifdef USEDOUBLE
#define real double
#else
#define real float
#endif

#ifdef I_WANT_TIME
#include <sys/times.h>
/*#include <limits.h>*/
#include <time.h>
long getRuntime(void) ;
#endif

/// Sparse definitions //////////////
typedef struct
{
    int indice;
    real valeur;
} sreal;

int sparseLineLength(sreal *line);
/////////////////////////////////////

#define VERSION "SVMTorch II V1.77 [Ennastie]"

#endif
