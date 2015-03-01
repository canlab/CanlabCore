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


#include "Cache.h"
#include "SVM.h"
#include "Kernel.h"

Cache::Cache(SVM &svm_, real taille_en_megs) : svm(svm_), kernel(*svm_.kernel)
{
  // Alloc
  l = svm.l;
  taille = (int)(taille_en_megs*1048576./((real)sizeof(real)*l));
  index_dans_liste = new Liste *[l];
  cached = new Liste[taille];
  cached_sauve = cached;

  cout << "# Max columns in cache: " << taille << "\n";
  if(taille < 2)
  {
    cout << "$ Change cache size : it's too small.\n\n";
    exit(0);
  }

  // Init
  Liste *ptr = cached;
  for(int i = 0; i < l; i++)
    index_dans_liste[i] = NULL;

  for(int i = 0; i < taille; i++)
  {
    ptr->adr = new real [l];
    ptr->index = -1;
    if(i != 0)
      ptr->prev = (ptr-1);
    else
      ptr->prev = &cached[taille-1];
    if(i != taille-1)
      ptr->suiv = (ptr+1);
    else
      ptr->suiv = cached;

    ptr++;
  }
}

void Cache::efface()
{
  Liste *ptr = cached;
  for(int i = 0; i < taille; i++)
  {
    ptr->index = -1;
    ptr = ptr->suiv;
  }

  for(int i = 0; i < l; i++)
    index_dans_liste[i] = NULL;
}

real *Cache::adresseCache(int index)
{
  Liste *ptr;

  // Rq: en regression faudrait faire gaffe a pas recalculer deux trucs...
  // Mais pb: -1 +1 a inverser dans la matrice...
  // Donc fuck.

  ptr = index_dans_liste[index];
  if( (ptr != NULL) && (ptr != cached) )
  {
//    cout << "Index " << index << " is already inside" << endl;
    ptr->prev->suiv = ptr->suiv;
    ptr->suiv->prev = ptr->prev;

    ptr->suiv = cached;
    ptr->prev = cached->prev;
    cached->prev->suiv = ptr;
    cached->prev = ptr;
    cached = ptr;
  }
  else
  {
    cached = cached->prev;
    if(cached->index != -1)
      index_dans_liste[cached->index] = NULL;
    cached->index = index;
    index_dans_liste[index] = cached;
    rempliColonne(index, cached->adr);
  }

//   {
//     cout << "Inside cache:" << endl;
//     Liste *ptx;
//     ptx = cached;
//     for(int i = 0; i < taille; i++)
//     {
//       cout << ptx->index << " ";
//       ptx = ptx->suiv;
//     }
//     cout << endl;
//     getchar();
//   }

  return(cached->adr);
}

Cache::~Cache()
{
  delete[] index_dans_liste;

  Liste *ptr = cached;
  for(int i = 0; i < taille; i++)
  {
    delete[] ptr->adr;
    ptr = ptr->suiv;
  }

  delete[] cached_sauve;  
}

CacheClassification::CacheClassification(SVM &svm_, real taille_en_megs)
  : Cache(svm_, taille_en_megs)
{
  y = svm.y;
}

void CacheClassification::rempliColonne(int index, real *adr)
{
  if(svm.deja_shrink && !svm.unshrink_mode)
  {
    if(y[index] > 0)
    {
      for(int it = 0; it < svm.n_active_var; it++)
      {
        int t = svm.active_var[it];
        adr[t] =  y[t]*kernel.evalue(index, t);
      }
    }
    else
    {
      for(int it = 0; it < svm.n_active_var; it++)
      {
        int t = svm.active_var[it];
        adr[t] = -y[t]*kernel.evalue(index, t);
      }
    }
  }
  else
  {
    if(y[index] > 0)
    {
      for(int i = 0; i < l; i++)
        adr[i] =  y[i]*kernel.evalue(index, i);
    }
    else
    {
      for(int i = 0; i < l; i++)
        adr[i] = -y[i]*kernel.evalue(index, i);
    }
  }
}

CacheRegression::CacheRegression(SVM &svm_, real taille_en_megs)
  : Cache(svm_, taille_en_megs)
{
  lm = l/2;
}

void CacheRegression::rempliColonne(int index, real *adr)
{
  int indexm = index%lm;
  if(svm.deja_shrink && !svm.unshrink_mode)
  {
    if(index < lm)
    {
      for(int i = 0; i < svm.n_active_var; i++)
      {
        int k = svm.active_var[i]%lm;
        adr[k] = kernel.evalue(indexm, k);
      }
    }
    else
    {
      for(int i = 0; i < svm.n_active_var; i++)
      {
        int k = svm.active_var[i]%lm;
        adr[k] = -kernel.evalue(indexm, k);
      }
    }
  }
  else
  {
    if(index < lm)
    {
      for(int i = 0; i < lm; i++)
        adr[i] =  kernel.evalue(indexm, i);
    }
    else
    {
      for(int i = 0; i < lm; i++)
        adr[i] = -kernel.evalue(indexm, i);
    }

  }

  for(int i = 0; i < lm; i++)
    adr[i+lm] = -adr[i];
}
