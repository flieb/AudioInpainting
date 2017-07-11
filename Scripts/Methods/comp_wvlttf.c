#include "mex.h"
#include "matrix.h"
#include "stdio.h"

/*
 * comp_wvltif.c - partial wavelet transform for speedup
 *
 * 
 *
 * December 2016, F. Lieb
 */
 

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  

  const mxArray *f = prhs[0];
  const double *g = mxGetPr(prhs[1]);
  const double *am = mxGetPr(prhs[2]);
  const mwSize *dims = mxGetDimensions(prhs[1]);
  mxArray *tmp;
  
  int Ls = (int) dims[0];
  int numwins = (int) dims[1];
  int m;
  
  double *tmp_r, *tmp_i;
  double *f_r = mxGetPr(f);
  double *f_i = mxGetPi(f);
  
  //printf("size: %d x %d\n",dims[0],dims[1]);
  
  //dims[0] is size(c,1), dims[1] is size(c,2)
  //printf("Number: %d\n", Ls);
  //printf("Number: %d\n", numwins);
 

  plhs[0] = mxCreateCellMatrix( (mwSize) numwins, (mwSize) 1);
  
  
  for (int i=0; i<numwins; ++i) {
      m = Ls/am[i];
      tmp = mxCreateDoubleMatrix((mwSize) m, (mwSize) 1, mxCOMPLEX );
      tmp_r = mxGetPr(tmp);
      tmp_i = mxGetPi(tmp);
      
      for(int j=0; j<am[i]; j++) {
          for (int k=0; k<m; ++k) {
              tmp_r[k] += (1/am[i])*f_r[ m*j + k ]*g[ (i*Ls)+(m*j + k) ];
              tmp_i[k] += (1/am[i])*f_i[ m*j + k ]*g[ (i*Ls)+(m*j + k) ];
          }
      }
      
      
      mxSetCell(plhs[0],i,mxDuplicateArray(tmp));
      
  }
  mxDestroyArray(tmp);
  
}
