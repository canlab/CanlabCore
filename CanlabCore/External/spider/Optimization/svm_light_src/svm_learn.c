
#include <mex.h> 
#include "svm_common.h"
char docfile[200];           /* file with training examples */
char modelfile[200];         /* file for resulting classifier */
/* interface to QP-solver */
double *optimize_qp(QP *, double *, long, double *, LEARN_PARM *);
int nx = 400;                /* Maximum size of QP-subproblems */
void   svm_learn(DOC *, long *, long, long, LEARN_PARM *, KERNEL_PARM *, 
		 KERNEL_CACHE *, MODEL *);
double compute_objective_function(double *, double *, long *, long *);
void   clear_index(long *);
void   add_to_index(long *, long);
long   compute_index(long *,long, long *);
void   optimize_svm(DOC *, long *, long *, long *, long *, MODEL *, long, 
		    long *, long, double *, double *, LEARN_PARM *, CFLOAT *, 
		    KERNEL_PARM *, QP *, double *);
void   compute_matrices_for_optimization(DOC *, long *, long *, long *, 
					 long *, long *, MODEL *, double *, 
					 double *, long, long, LEARN_PARM *, 
					 CFLOAT *, KERNEL_PARM *, QP *);
long   calculate_svm_model(DOC *, long *, long *, double *, double *, 
			   double *, LEARN_PARM *, long *, MODEL *);
long   check_optimality(MODEL *, long *, long *, double *, double *,long, 
			LEARN_PARM *,double *, double, long *, long *, long *,
			long *, long, KERNEL_PARM *);
long   identify_inconsistent(double *, long *, long *, long, LEARN_PARM *, 
			     long *, long *);
long   incorporate_unlabeled_examples(MODEL *, long *,long *, long *,
				      double *, double *, long, double *,
				      long *, long *, long, KERNEL_PARM *,
				      LEARN_PARM *);
void   update_linear_component(DOC *, long *, long *, double *, double *, 
			       long *, long, long, KERNEL_PARM *, 
			       KERNEL_CACHE *, double *,
			       CFLOAT *, double *);
long   select_next_qp_subproblem_grad(long *, long *, double *, double *, long,
				      long, LEARN_PARM *, long *, long *, 
				      long *, double *, long *, KERNEL_CACHE *,
				      long *, long *);
long   select_next_qp_subproblem_grad_cache(long *, long *, double *, double *,
					    long, long, LEARN_PARM *, long *, 
					    long *, long *, double *, long *,
					    KERNEL_CACHE *, long *, long *);
void   select_top_n(double *, long, long *, long);
long   shrink_problem(LEARN_PARM *, long *, long *,long *, long *, long *, 
		      long, long, long, double *, double **, long *);
void   reactivate_inactive_examples(long *, long *, double *, double **, 
				    double *, long, long, long, LEARN_PARM *, 
				    long *, long *, long *, long, DOC *, 
				    KERNEL_PARM *,KERNEL_CACHE *, MODEL *, 
				    CFLOAT *, double *, double *);
/* cache kernel evalutations to improve speed */
void   get_kernel_row(KERNEL_CACHE *,DOC *, long, long, long *, CFLOAT *, 
		      KERNEL_PARM *);
void   cache_kernel_row(KERNEL_CACHE *,DOC *, long, KERNEL_PARM *);
void   cache_multiple_kernel_rows(KERNEL_CACHE *,DOC *, long *, long, 
				  KERNEL_PARM *);
void   kernel_cache_shrink(KERNEL_CACHE *,long, long, long *);
void   kernel_cache_init(KERNEL_CACHE *,long, long);
void   kernel_cache_cleanup(KERNEL_CACHE *);
long   kernel_cache_malloc(KERNEL_CACHE *);
void   kernel_cache_free(KERNEL_CACHE *,long);
long   kernel_cache_free_lru(KERNEL_CACHE *);
CFLOAT *kernel_cache_clean_and_malloc(KERNEL_CACHE *,long);
long   kernel_cache_touch(KERNEL_CACHE *,long);
long   kernel_cache_check(KERNEL_CACHE *,long);
double estimate_margin_vcdim(MODEL *, double, double, KERNEL_PARM *);
double estimate_sphere(MODEL *, KERNEL_PARM *);
void   write_model(char *, MODEL *, KERNEL_PARM *);
void   write_prediction(char *, MODEL *, double *, double *, long *, long *,
			long, LEARN_PARM *);
void   write_alphas(char *, double *, long *, long);
void   read_input_parameters(int, char **, char *, char *,long *, long *, 
			     LEARN_PARM *, KERNEL_PARM *);
void   print_help();
typedef struct cache_parm_s {
  KERNEL_CACHE *kernel_cache;
  CFLOAT *cache;
  DOC *docs; 
  long m;
  KERNEL_PARM *kernel_parm;
  long offset,stepsize;
} cache_parm_t;
/*---------------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nhrs, const mxArray *prhs[])
  /*
int main (argc, argv1)
int argc;
char *argv1[]; */
{ 
  DOC *docs;  /* training examples */
  long *label,max_docs,max_words_doc;
  long totwords,totdoc,ll;
  long kernel_cache_size;
  KERNEL_CACHE kernel_cache;
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  MODEL model;
  double *R; int i,j,l,n;
  char *k,*argv[100] ={"svm","train.dat"};
  int argc=2;
  double *X,*Y,*C;
  double np,nt;
  char   ker; double p1;
  double nrm;

  read_input_parameters(argc,argv,docfile,modelfile,&verbosity,
			&kernel_cache_size,&learn_parm,&kernel_parm);
  /* READ MEX */
  X = (double *) mxGetPr(prhs[0]);
  Y = (double *) mxGetPr(prhs[1]);
  l = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  C = (double *) mxGetPr(prhs[2]);
  learn_parm.svm_c=C[0];
  learn_parm.biased_hyperplane=1;
  printf("l=%d  n=%d\n",l,n);
  // kernel_parm.mycustom: what they ,mean
  // 0=type of kernel,  1=ridge,2=bla ridge const,  3,4=bal ridge n+/n,n-/n
  
  C = (double *) mxGetPr(prhs[3]);kernel_parm.mycustom[1]=C[0];//ridge
  C = (double *) mxGetPr(prhs[4]);kernel_parm.mycustom[2]=C[0];//bal ridge
  C = (double *) mxGetPr(prhs[5]);i=(int) C[0];
  kernel_parm.mycustom[0]=i; // type of kernel
  //if(i==1) ker='l'; if(i==2) ker='p'; if(i==3) ker='r'; 
  ker='s'; // special kernel implements all kernels - lin,poly,rbf
  C = (double *) mxGetPr(prhs[6]);p1=C[0];
  //printf("%f %c %f\n",kernel_parm.mycustom[3],ker,p1);
  if(ker=='l') kernel_parm.kernel_type=0;
  if(ker=='p') kernel_parm.kernel_type=1;
  if(ker=='r') kernel_parm.kernel_type=2;
  if(ker=='s') kernel_parm.kernel_type=4;
  //type of kernel
  if(i==1) //linear
    {
      kernel_parm.coef_lin=1.0; //1.0/n;
      kernel_parm.coef_const=0; 
      p1=1; // linear, no powers
  }
  else//poly or other
    {
      kernel_parm.coef_lin=1.0; //1.0/n;
      kernel_parm.coef_const=1; 
    }
  np=0; for(i=0;i<l;i++) if(Y[i]==1) np++;
  nt=((double)l);
  kernel_parm.mycustom[3]=np/nt;
  kernel_parm.mycustom[4]=(nt-np)/nt;
  kernel_parm.poly_degree = p1;
  kernel_parm.rbf_gamma = p1;
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
  label = (long *)my_malloc(sizeof(long)*(l+10));  
  for(i=0;i<l;i++) 
    {
       label[i]=Y[i];   
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
       docs[i].label=Y[i]; 
       docs[i].words[n].wnum=0; 
       docs[i].words[n].weight=0; 
    }
  totwords=n;
  totdoc=l;
 
  if(kernel_parm.kernel_type == LINEAR) { /* don't need the cache */
    svm_learn(docs,label,totdoc,totwords,&learn_parm,&kernel_parm,
	      NULL,&model);
  }
  else {
    kernel_cache_init(&kernel_cache,totdoc,kernel_cache_size);
    svm_learn(docs,label,totdoc,totwords,&learn_parm,&kernel_parm,
	      &kernel_cache,&model);
    kernel_cache_cleanup(&kernel_cache);
  }
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
  // indexes
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
  free(label);
}


