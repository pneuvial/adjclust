//C functions to make diagonal elements 1 by using transformation s(i,j) = s(i,j) + 1 - 0.5*(s(i,i) + s(j,j))

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP CMakeDiagOne(SEXP X){
  
  int p, i, j;
  double *xptr,*outptr;
  xptr = REAL(X);
  p = INTEGER(GET_DIM(X))[0];

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, p, p));
  outptr = REAL(out);

  for ( i = 0; i < p; i++ ) {
    for ( j = 0; j < p ; j++ ) {
      outptr[i + j*p] = xptr[i + j*p] + 1 - 0.5*(xptr[i + i*p] + xptr[j + j*p]);
    }
  }

  UNPROTECT(1);
  return(out);
}
