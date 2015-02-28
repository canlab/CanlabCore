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


#ifndef INC_CACHE
#define INC_CACHE

#include "general.h"
//#include "SVM.h"
//#include "Kernel.h"

class SVM;
class Kernel;

typedef struct Liste_
{
    real *adr;
    int index;
    Liste_ *prev;
    Liste_ *suiv;
} Liste;

class Cache
{
  private:
  int taille;
  Liste *cached, *cached_sauve;
  Liste **index_dans_liste;
  
  protected:
  int l;
  SVM &svm;
  Kernel &kernel;

  virtual void rempliColonne(int index, real *adr) = 0;

  public:
  Cache(SVM &svm_, real taille_en_megs);
  real *adresseCache(int index);  
  void efface();

  virtual ~Cache();
};

class CacheClassification : public Cache
{
  private:
  real *y;

  public:
  CacheClassification(SVM &svm_, real taille_en_megs);
  virtual void rempliColonne(int index, real *adr);
};

class CacheRegression : public Cache
{
  private:
  int lm;

  public:
  CacheRegression(SVM &svm_, real taille_en_megs);
  virtual void rempliColonne(int index, real *adr);
};

#endif
