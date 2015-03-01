/*

// A MEX interface to SVM Light - call with 
// [a b0 xsv]=svm_light(X,Y,C,ridge,balr,ker_type,ker_val,verb,[argv],[len])
//
// X is a  examples x features matrix
// Y is a example outputs x 1 matrix 
// C is the C soft margin parameter
// ridge is the value of ridge on the kernel matrix
// balr is another version of the ridge balanced by the 
//   number of pos and neg examples [NOTE: this isn't implemented yet]
// ker_type=1,2,3 or 4 for 1:linear,2:polynomial,3:rbf or 4:special (custom)
//   if custom, the X data should be the kernel matrix instead of the data
// ker_val is the value of rbf sigma or the polynomial degree
// verb is 0 for no text from svmlight or >0 for text
// these last two options can be omitted:
// argv -a string (can't be an empty string if parameter exists) 
//       of options e.g '-i1 -v2' to pass to svmlight
//       if you want no options use argv=' '
// len - length(argv)
//
// [a] are the alphas for vectors with indexes [xsv]
// and b0 is the threshold
//

// e.g a typical matlab call:
// X=rand(100,5); Y=sign(sum(X,2)-2.5); %% make data
// [a b c]= svm_light(X,Y,1000,0,0,'l',1,1);  %% linear kernel on problem
*/

#include <mex.h> 
# include "svm_common.h"
# include "svm_learn.h"
char docfile[200];           /* file with training examples */
char modelfile[200];         /* file for resulting classifier */




double mycustom[10]; /* used for storing ridge, etc. in custom kernel*/



void mexFunction(int nlhs, mxArray *plhs[], int nhrs, const mxArray *prhs[])
{ 
  DOC *docs;  /* training examples */
  long max_docs,max_words_doc;
  long totwords,totdoc,ll;
  long kernel_cache_size;
  double *target;
  KERNEL_CACHE kernel_cache;
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  MODEL model;
  /* my stuff below */
  double *R; int i,j,l,n;
  char *k,*k2,*argv[100] ={"svm",""};
  int argc=4;
  double *X,*Y,*C;
  double np,nt;
  char   ker; double p1;
  double nrm;

  n=0;
  if (nhrs>8)
  {
    C = (double *) mxGetPr(prhs[9]);i=(int) C[0]+1; /* argv string size*/
    k=malloc(sizeof(char)*i+1); 
    mxGetNChars(prhs[8], k, i); k[i]='\0';
    l=0; 
    k2=malloc(sizeof(char)*i+1); 
    for (j=0;j<i;j++)
    {
      k2[l]=k[j];
      if ( l>0 & (k[j]==' ' | j==i-1)) /* new argv*/
	{
	  k2[l]='\0'; argv[1+n]=k2;
	  l=0;n++;  k2=malloc(sizeof(char)*i+1); 
	}
      else { l++; }
    }
  }
  k2=malloc(sizeof(char)*10+1); k2[0]='j'; k2[1]='\0'; 
  argv[n+1]=k2; 
  k2=malloc(sizeof(char)*10+1); k2[0]='j'; k2[1]='\0'; 
  argv[n+2]=k2;
  argc=3+n;
    

  read_input_parameters(argc,argv,docfile,modelfile,&verbosity,
			&kernel_cache_size,&learn_parm,&kernel_parm);
 
  for(i=1;i<argc;i++)
     free(argv[i]);
  if (nhrs>8)
    {
      free(k);
    }

  /* READ MEX */
  X = (double *) mxGetPr(prhs[0]);
  Y = (double *) mxGetPr(prhs[1]);
  l = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  C = (double *) mxGetPr(prhs[2]);
  learn_parm.svm_c=C[0];;
  learn_parm.biased_hyperplane=1;

  /* kernel_parm.custom: what they ,mean
   0=type of kernel,  1=ridge,2=bla ridge const,  3,4=bal ridge n+/n,n-/n
  */

  C = (double *) mxGetPr(prhs[3]);
  mycustom[1]=C[0];/*ridge*/
  C = (double *) mxGetPr(prhs[4]);
  mycustom[2]=C[0];/*bal ridge*/
  C = (double *) mxGetPr(prhs[5]);i=(int) C[0];
  mycustom[0]=i; /* type of kernel */
   ker='s'; /* special kernel implements all kernels - lin,poly,rbf */
  C = (double *) mxGetPr(prhs[6]);p1=C[0];
  C = (double *) mxGetPr(prhs[7]);verbosity=C[0];

  if(ker=='l') { kernel_parm.kernel_type=0; }
  if(ker=='p') { kernel_parm.kernel_type=1; }
  if(ker=='r') { kernel_parm.kernel_type=2; }
  if(ker=='s') { kernel_parm.kernel_type=4; }
  /*type of kernel */
  if(i==1) /*linear */
    {
      kernel_parm.coef_lin=1.0; /*1.0*/;
      kernel_parm.coef_const=0; 
      p1=1; /* linear, no powers*/
  }
  else/*poly or other*/
    {
      kernel_parm.coef_lin=1.0; /*1.0 */
      kernel_parm.coef_const=1; 
    }
  np=0; for(i=0;i<l;i++) if(Y[i]==1) np++;
  nt=((double)l);
  mycustom[3]=np/nt;
  mycustom[4]=(nt-np)/nt;
  kernel_parm.poly_degree = p1;
  kernel_parm.rbf_gamma = p1;

  if (verbosity>0 ) 
    printf("l=%d  n=%d\n",l,n);
  /*
  if(verbosity>=1) {
    printf("Scanning examples..."); fflush(stdout);
  }
  nol_ll(docfile,&max_docs,&max_words_doc,&ll); 
  max_words_doc+=2;
  ll+=2;
  max_docs+=2;
  if(verbosity>=1) {
    printf("done\n"); fflush(stdout);
  } 
  */
  /*
  read_documents(docfile,docs,label,max_words_doc,ll,&totwords,&totdoc);
  */
 
  docs = (DOC *)my_malloc(sizeof(DOC)*(l+10));
  target = (double *)my_malloc(sizeof(double)*(l+10));  
  for(i=0;i<l;i++) 
    {
       target[i]=Y[i];   
       docs[i].words = (WORD *)my_malloc(sizeof(WORD)*(n+1));
       docs[i].docnum=i;
       
       nrm=0;
       for(j=0;j<n;j++) 
	 {
	   docs[i].words[j].wnum=j+1; 
	   docs[i].words[j].weight=X[i+j*l]; 
	   nrm+=X[i+j*l]*X[i+j*l]; 
	 } 
       docs[i].twonorm_sq=nrm;
       docs[i].words[n].wnum=0; 
       docs[i].words[n].weight=0; 
    }
  totwords=n;
  totdoc=l;
 
  if(kernel_parm.kernel_type == LINEAR) { /* don't need the cache */
    if(learn_parm.type == CLASSIFICATION) {
      svm_learn_classification(docs,target,totdoc,totwords,&learn_parm,
			       &kernel_parm,NULL,&model);
    }
    else if(learn_parm.type == REGRESSION) {
      svm_learn_regression(docs,target,totdoc,totwords,&learn_parm,
			   &kernel_parm,NULL,&model);
    }
  }
  else {
    if(learn_parm.type == CLASSIFICATION) {
      /* Always get a new kernel cache. It is not possible to use the
         same cache for two different training runs */
      kernel_cache_init(&kernel_cache,totdoc,kernel_cache_size);
      svm_learn_classification(docs,target,totdoc,totwords,&learn_parm,
			       &kernel_parm,&kernel_cache,&model);
      /* Free the memory used for the cache. */
          kernel_cache_cleanup(&kernel_cache);
    }
    else if(learn_parm.type == REGRESSION) {
      /* Always get a new kernel cache. It is not possible to use the
         same cache for two different training runs */
      kernel_cache_init(&kernel_cache,2*totdoc,kernel_cache_size);
      svm_learn_regression(docs,target,totdoc,totwords,&learn_parm,
			   &kernel_parm,&kernel_cache,&model);
      /* Free the memory used for the cache. */
      kernel_cache_cleanup(&kernel_cache);
    }
  }

  /* -------------------------------

  /*   write_model(modelfile,&model,&kernel_parm); */
  /* WRITE MEX */
  /*  printf("number of poos : %d %d\n",l,model.sv_num); */
  /*
  plhs[0] = mxCreateDoubleMatrix(l,1,mxREAL);
  R=mxGetPr(plhs[0]);
  for(i=0;i<l;i++) R[i]=0;
  for(i=0;i<l;i++) 
    if(model.index[i]>=1 && model.index[i]<=l)
    R[i]=model.alpha[model.index[i]];
  */
  plhs[0] = mxCreateDoubleMatrix(model.sv_num-1,1,mxREAL);  
  R=mxGetPr(plhs[0]);
  for(i=0;i<model.sv_num-1;i++) 
    R[i]=model.alpha[i+1];
  
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  R=mxGetPr(plhs[1]);
  R[0]=-model.b;
  
  /*
  plhs[2] = mxCreateDoubleMatrix(model.sv_num-1,n,mxREAL);
  R=mxGetPr(plhs[2]);
  for(j=1;j<(model.sv_num);j++)
   for(i=0;i<n;i++) 
     {
      R[(i)*((model.sv_num)-1)+(j-1)]= model.supvec[j]->words[i].weight; 
     }
*/
  /* indexes */
  plhs[2] = mxCreateDoubleMatrix(model.sv_num-1,1,mxREAL);  
  R=mxGetPr(plhs[2]);
  for(i=0;i<model.sv_num-1;i++) 
      R[i]=(double) model.supvec[i+1]->docnum;
  /* WRITE MEX DONE */
  free(model.supvec);
  free(model.alpha);
  free(model.index);
  for(i=0;i<l;i++) 
      free(docs[i].words);
  
  free(docs);
  free(target);
}


/* ------------------------*/


double get_x(a,i) /* get index i in example a , used for custom kernel*/
WORD *a; int i;          
{
    register FVAL sum=0;
    register WORD *ai;
    ai=a;
    while(ai->wnum < i) ai++;
    sum=ai->weight;
    return((double)sum);
}
