//C functions

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>


int min(int a,int b)
{
  if (a <= b)  return a;
  else  return b; 
} 

SEXP CSparseMatL(SEXP X, SEXP Rp, SEXP Rh)
{
  int *p, *h, i ,j;
  double *xptr, *outptr;
  xptr = REAL(X);
  p = INTEGER(Rp);
  h = INTEGER(Rh);

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, *p, *h));
  outptr = REAL(out);

  // intialize "out" matrix values to zero
  for(i = 0; i < *p; i++)
  {
    for(j = 0; j < *h ; j++)
    {
        outptr[i + (*p)*j] = 0;      
    }
  }

  for(i = 0; i < *p; i++)
  {
    for(j = 0; j < *h ; j++)
    {
      if( (j < *h)&&(i+j+1 < *p) )
        outptr[i + (*p)*j] = xptr[i+j+1 +(*p)*(j+1)];      
    }
  }

  UNPROTECT(1);
  return(out);
}


SEXP CSparseRmatR(SEXP X, SEXP Rp, SEXP Rh)
{
  int *p, *h, i ,j;
  double *xptr, *outptr;
  xptr = REAL(X);
  p = INTEGER(Rp);
  h = INTEGER(Rh);

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, *p, *h) );
  outptr = REAL(out);

/*
  // intialize "out" matrix values to zero
  for(j = 0; j < *p; j++)
  {
    for(k = 0; k < *h ; k++)
    {
      outptr[j + (*p)*k] = 0;      
    }
  }

	l=0;
  for(i = 1; i < *p ; i++)
  {
    j = i;
    k = 1;
    while( (j<*p)&&(k<*h+1) )
    {
      tempptr[l++] = xptr[j + (*p)*k];
      j++;k++;
    }
  }
  l--;

  for(i = 0; i < *p ; i++)
  {
  	for(j = 0 ;j < min(*h,*p-1-i); j++)
  	{
      outptr[i + (*p)*j] = tempptr[l--];
  	}
  }
*/
  
  for(i = *p-1 ; i >= 0 ; i--)
  {
    for(j = 1; j <= *h ; j++)
    {
        outptr[(*p - 1 - i) + (*p)*(j-1)] = xptr[i + (*p)*(j)];      
    }
  }


  UNPROTECT(1);
  return(out);

}
