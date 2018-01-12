//C function to check condition s(i,j) <= 0.5*(s(i,i) + s(j,j)) in band derived from sparse matrix

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP CSparseCondnCheck(SEXP X, SEXP Rp, SEXP Rh){
  int i, j, *p, *h;  
  double *xptr;
  xptr = REAL(X);
  p = INTEGER(Rp);
  h = INTEGER(Rh);

  SEXP out = PROTECT(out = allocVector(LGLSXP, 1));

  LOGICAL(out)[0] = TRUE;

  for (i = 0; i < *p; i++) {
    for (j = 0; (j < *h+1)&&(j < i+1); j++) {
      if (xptr[i + (*p)*j] > 0.5*(xptr[i] + xptr[i-j])){
        LOGICAL(out)[0] = FALSE;
        i = *p + 1; //to break outer loop
        break;
      }
    }
  }
  
  UNPROTECT(1);
  return(out);
}