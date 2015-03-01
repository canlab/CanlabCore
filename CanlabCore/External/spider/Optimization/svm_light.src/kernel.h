/************************************************************************/
/*                                                                      */
/*   kernel.h                                                           */
/*                                                                      */
/*   User defined kernel function. Feel free to plug in your own.       */
/*                                                                      */
/*   Copyright: Thorsten Joachims                                       */
/*   Date: 16.12.97                                                     */
/*                                                                      */
/************************************************************************/
/* KERNEL_PARM is defined in svm_common.h The field 'mycustom' is reserved for */
/* parameters of the user defined kernel. You can also access and use */
/* the parameters of the other kernels. */

extern double mycustom[10];
extern double get_x(WORD *a, int i);


double custom_kernel(kernel_parm,a,b) /* plug in you favorite kernel */
KERNEL_PARM *kernel_parm;      
DOC *a,*b;
{
  double ridge=0;
  if(a==b)
    {
      /* if(a->label==1)  // unbalanced ridge not implemented yet
	ridge= mycustom[3];
      else 
	ridge= mycustom[4];
	ridge *= mycustom[2]; // balanced ridge */
      ridge=ridge+mycustom[1];    
    }

  
  if ((mycustom[0])<3) /* polynomial */
    {
    return ridge+((CFLOAT)pow(kernel_parm->coef_lin*sprod_ss(a->words,b->words)+kernel_parm->coef_const,(double)kernel_parm->poly_degree)); 
    }
  
  if ((mycustom[0])==3) /* RBF */
    {
  return  ridge+((CFLOAT)exp(-(a->twonorm_sq-2*sprod_ss(a->words,b->words)+b->twonorm_sq)/kernel_parm->rbf_gamma));
    }

  if ((mycustom[0])==4) /* custom */ 
      {
	return  ridge+get_x(b->words,1+get_x(a->words,1)) ;
      }
}









