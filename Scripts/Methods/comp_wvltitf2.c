#include "mex.h"
#include "stdio.h"

/*
 * comp_wvltif.c - partial inverse wavelet transform for speedup
 *
 * 
 *
 * December 2016, F. Lieb
 */
 

/*void timestwo(double y[], double x[])
{
  y[0] = 2.0*x[0];
}*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  

  const mxArray *c = prhs[0];
  const double *gd = mxGetPr(prhs[1]);
  const double *am = mxGetPr(prhs[2]);
  const mwSize *dims = mxGetDimensions(prhs[1]);
  //const mxArray *cellArray = mxGetCell(c,0);
  //const mwSize *dim2 = mxGetDimensions(cellArray);
  double *data_r,*data_i,*data_g,*zr,*zi;
  int Ls = (int) dims[0];
  int numwins = (int) dims[1];
  int indx = 0;
  size_t m;
  
  //dims[0] is size(c,1), dims[1] is size(c,2)
  //printf("Number: %d\n", Ls);
  //printf("Number: %d\n", numwins);
 
  
  plhs[0] = mxCreateDoubleMatrix( (mwSize) 1, (mwSize) Ls, mxCOMPLEX);
  zr = mxGetPr(plhs[0]);
  zi = mxGetPi(plhs[0]);
  
  //dims[0]
  for(int i=0; i<numwins; ++i) {
      data_r = mxGetPr(mxGetCell(c,i));
      data_i = mxGetPi(mxGetCell(c,i));
      m = mxGetM(mxGetCell(c,i));
      //printf("length: %d",m);
      for (int j=0; j<Ls; ++j) {
          //indx = j % m;
          zr[j] += gd[ (i*Ls)+j ]*data_r[j%m];
          zi[j] += gd[ (i*Ls)+j ]*data_i[j%m];
          //zr[j] += gd[ (j*numwins)+i ]*data_r[j%m];
          //zi[j] += gd[ (j*numwins)+i ]*data_i[j%m];
      }
  }
  
  
  /* Create matrix for the return argument. */
  //plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
  
  /* Assign pointers to each input and output. */
  //x = mxGetPr(prhs[0]);
  //y = mxGetPr(plhs[0]);
  
  /* Call the timestwo subroutine. */
  //timestwo(y,x);
}
