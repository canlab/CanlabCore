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

int sparseLineLength(sreal *line)
{
  int i = 0;
  while(line->indice != -1)
  {
    i++;
    line++;
  }
  return(i);
}

#ifdef I_WANT_TIME
long getRuntime(void) 
{
        struct tms buffer;
        times(&buffer);
        return((long)(buffer.tms_utime+buffer.tms_stime));
}
#endif
