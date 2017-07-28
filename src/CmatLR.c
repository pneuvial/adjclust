//C functions to create matL and RmatR matrix

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>


int minm(int a,int b)
{
  if (a <= b)  return a;
  else  return b; 
} 

SEXP CmatL(SEXP X, SEXP Rp, SEXP Rh)
{
  int *p, *h, i ,j, k;
  double *xptr, *outptr;
  xptr = REAL(X);
  p = INTEGER(Rp);
  h = INTEGER(Rh);

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, *p, *h));
  outptr = REAL(out);

  // intialize "out" matrix values to zero
  for(j = 0; j < *p; j++)
  {
    for(k = 0; k < *h ; k++)
    {
      outptr[j + (*p)*k] = 0;      
    }
  }

  for(i=0; i < *p; i++)
  {
  	for(j = i+1 ; j <= minm(i + *h,*p - 1); j++)
  	{
      outptr[i + (*p)*(j-i-1)] = xptr[i + (*p)*j];
  	}
  }

  UNPROTECT(1);
  return(out);
}


SEXP CRmatR(SEXP X, SEXP Rp, SEXP Rh)
{
  int *p, *h, j, k;
  double *xptr, *outptr;
  xptr = REAL(X);
  p = INTEGER(Rp);
  h = INTEGER(Rh);

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, *p, *h) );
  outptr = REAL(out);

  // intialize "out" matrix values to zero
  for(j = 0; j < *p; j++)
  {
    for(k = 0; k < *h ; k++)
    {
      outptr[j + (*p)*k] = 0;      
    }
  }

  for(j = 0; j < *p; j++)
  {
    for(k = 0; k < *h ; k++)
    {
      if( (*p-2-k-j) >= 0 )
      { 
      outptr[j + (*p)*k] = xptr[(*p-1-j) + (*p)*(*p-2-k-j)];
      }      
    }
  }
  
  
  UNPROTECT(1);
  return(out);

}
