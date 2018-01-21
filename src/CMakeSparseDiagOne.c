//C function to make diag one in band derived from sparse matrix

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP CMakeSparseDiagOne(SEXP X, SEXP Rp, SEXP Rh) {
  int i, j, *p, *h;  
  double *xptr, *outptr;
  xptr = REAL(X);
  p = INTEGER(Rp);
  h = INTEGER(Rh);

  SEXP out = PROTECT(out = allocMatrix(REALSXP, *p, *h+1));

  outptr = REAL(out);

  //initialize to zero  
  for(i = 0; i < *p; i++) {
    for(j = 0; j < *h+1; j++) {
      outptr[i + (*p)*j] =  0;
    }
  }

  for(i = 0; i < *p; i++) {
    for(j = 0; (j < *h+1)&&(j < i+1); j++) {
      outptr[i + (*p)*j] =  1 + xptr[i + (*p)*j] - 0.5*( xptr[i] + xptr[i-j] );
    }
  }
  
  UNPROTECT(1);
  return(out);
}
