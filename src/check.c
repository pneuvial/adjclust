//C function to check condition s(i,j) <= 0.5*(s(i,i) + s(j,j))

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP CCondnCheck(SEXP X){
  int p,i,j;  
  double *xptr;
  xptr = REAL(X);
  p = INTEGER(GET_DIM(X))[0];

  SEXP out = PROTECT(out = allocVector(LGLSXP, 1));

  LOGICAL(out)[0] = TRUE;

  for ( i = 0; i < p; i++ ) {
    for ( j = 0; j <= i ; j++ ) {
      if (xptr[i + j*p] > 0.5*(xptr[i + i*p] + xptr[j + j*p])) {
        LOGICAL(out)[0] = FALSE;
        j = p + 1; //to break outer for loop as well
        break; //breaks inner for loop
      }
    }
  }

  UNPROTECT(1);
  return(out);
}

