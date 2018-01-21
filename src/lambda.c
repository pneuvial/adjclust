//C function to find lambda diven a sparse band

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP CFindLambda(SEXP X, SEXP Rp, SEXP Rh){
  int i, j, *p, *h;  
  double *xptr, *tempptr;
  xptr = REAL(X);
  p = INTEGER(Rp);
  h = INTEGER(Rh);

  SEXP out = PROTECT(out = allocVector(REALSXP, 1));
  SEXP temp = PROTECT(temp = allocMatrix(REALSXP ,*p ,*h+1));
  tempptr = REAL(temp);
  REAL(out)[0] = 0;

  for(i = 0; i < *p; i++) {
    for(j = 0; (j < *h+1)&&(j < i+1); j++) {
      tempptr[i + (*p)*j] = xptr[i + (*p)*j] - xptr[i-j];
    }
  }

//  finding maximum
  for(i = 0; i < *p; i++) {
    for(j = 0; (j < *h+1)&&(j < i+1); j++) {
      if(tempptr[i + (*p)*j] > REAL(out)[0]) REAL(out)[0] = tempptr[i + (*p)*j];
    }
  }

  UNPROTECT(2);
  return(out);
}
