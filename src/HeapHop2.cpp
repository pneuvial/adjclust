#include <Rcpp.h>

#include "ClassesHeap.h"
#include "FusionnedClasses.h"
#include "PseudoMatrix.h"
#include <fstream>

// [[Rcpp::export]]
Rcpp::NumericVector  HeapHop2(Rcpp::NumericVector Input,const int p, const int h, const int NbCLasses){
  Rcpp::NumericVector merge=Rcpp::NumericVector(Rcpp::Dimension(3,p-1));
  // Rcpp::CharacterVector output(1);
  double *W = Input.begin();
  ClassesHeap H(W, int(p), int(h), int(NbCLasses));
  //  std::cout<<H;
  for(int index=0; index< 3*(p-1); index++) merge[index]=H.Output[index];
  return merge;
}

