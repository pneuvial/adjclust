//C functions to extract data around diagonal

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP DiagBand(SEXP X, SEXP h){
  int *hptr,p,len,i,j,k;
  double *xptr,*outptr;
  xptr = REAL(X);
  hptr = INTEGER(h);
  p = INTEGER(GET_DIM(X))[0];
  len = (p-1)*(*hptr) - ((*hptr)*(*hptr - 1)/2);

  SEXP out;
  PROTECT(out = allocVector(REALSXP,len));
  outptr = REAL(out);

  k=0;
  for ( i = 0; i < p; i++ ) {
    for ( j = ( i-(*hptr) > 0 ? i-(*hptr): 0 ); j< i; j++ ) {      
      outptr[k++] = xptr[i+j*p];
    }
  }

  UNPROTECT(1);
  return(out);
}

SEXP OneDiagBand(SEXP X, SEXP h){
  int *hptr,p,len,i,j,k;
  double *xptr,*outptr;
  xptr = REAL(X);
  hptr = INTEGER(h);
  p = INTEGER(GET_DIM(X))[0];
  len = (p-1)*(*hptr) - ((*hptr)*(*hptr - 1)/2) + p;

  SEXP out;
  PROTECT(out = allocVector(REALSXP,len));
  outptr = REAL(out);

  k=0;
  for ( i = 0; i < p; i++ ) {
    for ( j = i; j< i+(*hptr)+1 && j<p ; j++ ) { //Additonal one's included
      if (j == i) {
      outptr[k++] = 1;
      } else {
      outptr[k++] = xptr[i+j*p];
      }
    }
  }

  UNPROTECT(1);
  return(out);
}
