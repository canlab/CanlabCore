#include <mexarg.h>
#include <stdlib.h>
#include <string.h>

         int mexarg::isCellmatrix(int nrhs, const mxArray *prhs[])
         {
            return nrhs==1;
         }
    
        int mexarg::require(int nrhs, const mxArray *prhs[], const char *name)
        {
             if(! isCellmatrix(nrhs,prhs) ) 
             {
                for(int k=0;k<nrhs;k++)
                {
                    if(mxIsCell(prhs[k]))
                    {
                        char Buffer[256];
                        mxArray * cname;
                        cname= mxGetCell(prhs[k],0);
                        mxGetString(cname,Buffer,256);
                        if(strcmp(Buffer,name)==0)
                            return 1;
                    }
                }
             }
             return 0;
        }
        
        int mexarg::getstring(int nrhs, const mxArray *prhs[], const char *name,char *target,int len)
        {
            if(! isCellmatrix(nrhs,prhs) ) 
            {
                for(int k=0;k<nrhs;k++)
                {
                    if(mxIsCell(prhs[k]))
                    {
                        char Buffer[256];
                        mxArray * cname;
                        cname= mxGetCell(prhs[k],0);
                        mxGetString(cname,Buffer,256);
                        if(strcmp(Buffer,name)==0)
                        {
                            cname= mxGetCell(prhs[k],1);
                            mxGetString(cname,target,len);
                            return 1;
                        }
                    }
                }
            }
            return 0;
        }
        int mexarg::getmatrix(int nrhs,const mxArray *prhs[], const char *name,mxArray **target)
        {
            if(! isCellmatrix(nrhs,prhs) ) 
            {
                   for(int k=0;k<nrhs;k++)
                {
                    if(mxIsCell(prhs[k]))
                    {
                        char Buffer[256];
                        mxArray * cname;
                        cname= mxGetCell(prhs[k],0);
                        mxGetString(cname,Buffer,256);
                        if(strcmp(Buffer,name)==0)
                        {
                            *target= mxGetCell(prhs[k],1);
                            return 1;
                        }
                    }
                }
                   
            }
            
            return 0;
        }
        int mexarg::getscalar(int nrhs,const mxArray *prhs[], const char *name,double *target)
        {
            if( !isCellmatrix(nrhs,prhs) ) 
            {
                
                for(int k=0;k<nrhs;k++)
                {
                    if(mxIsCell(prhs[k]))
                    {
                        char Buffer[256];
                        mxArray * cname;
                        cname= mxGetCell(prhs[k],0);
                        mxGetString(cname,Buffer,256);
                        if(strcmp(Buffer,name)==0)
                        {
                            *target= mxGetScalar(mxGetCell(prhs[k],1));
                            return 1;
                        }
                    }
                }
            }
            return 0;
        }