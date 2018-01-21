//C functions to extract diagonal and data around it in upper matrix

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>


SEXP CSparseBand(SEXP Rxslot, SEXP Rpslot, SEXP Rislot, SEXP Rp, SEXP Rh){
  int *h, *p, l, q, count, arrind, j, k, *pslot, *islot;
  double *xslot, *outptr;
  xslot = REAL(Rxslot);
  pslot = INTEGER(Rpslot);
  islot = INTEGER(Rislot);
  h = INTEGER(Rh);
  p = INTEGER(Rp);

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, *p, *h+1));
  outptr = REAL(out);

  // intialize "out" matrix values to zero
  for(j = 0; j < *p; j++)
  {
    for(k = 0; k < *h+1 ; k++)
    {
      outptr[j + (*p)*k] = 0;      
    }
  }

  l = 0;
  for(j = 0; j < *p; j++)
  {
    count = pslot[j+1] - pslot[j];
    if(count == 0) continue;

    for(q = l; q <= l+count-1 ; q++)
    {
      arrind = j - islot[q] ;
      if( (arrind >= 0) && (arrind < *h+1) )
      {
        outptr[j + (*p)*arrind] = xslot[q];
      }     
    }
    l += count;
  }

  UNPROTECT(1);
  return(out);
}
