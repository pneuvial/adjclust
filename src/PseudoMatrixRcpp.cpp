#include <Rcpp.h>

#include "ClassesHeap.h"
#include "FusionnedClasses.h"
#include "PseudoMatrix.h"
#include <fstream>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector  HeapHop(Rcpp::NumericVector Input,const int p, const int h, const int NbCLasses){
  Rcpp::NumericVector merge=Rcpp::NumericVector(Rcpp::Dimension(3,p-1));
  // Rcpp::CharacterVector output(1);
  double *W = Input.begin();
  double *V = new double[p * (h+1)];
  //for (int i = 0;  i < p * (h + 1); i++)
  //  V[i] = 0;
  for (int Ligne = 0; Ligne < p - h; Ligne++)
    for (int Colonne = 0; Colonne < h+1; Colonne++)
      V[Ligne * (h + 1) + Colonne] = W[Ligne * (h + 1) + Colonne];
  for (int Ligne = p - h ; Ligne < p; Ligne++)
    for (int Colonne = 0; Colonne < h + 1; Colonne++)
      V[Ligne * (h+1) + Colonne] = 0;
  int Indice = (h+1) * (p - h);
  for (int Ligne = p - h ; Ligne < p; Ligne++)
    for (int Colonne = 0; Colonne < h - (Ligne - p + h); Colonne++)
      V[Ligne * (h + 1) + Colonne] = W[Indice++];

/*  
  for (int Ligne = 0; Ligne < p; Ligne++) {
    for (int Colonne = 0; Colonne < h+1; Colonne++) {
      Rcpp::Rcout << V[Ligne  * (h+1) + Colonne] << '\t';
    }
    Rcpp::Rcout << std::endl;
  }
*/
  // for (int Index = 0; Index<  p * (h+1) - (h*(h+1)/2); Index++) {
  //   Rcpp::cout << W[Index] << '\t';
  //   Rcpp::cout << std::endl;
  // }


  // std::ofstream OutputFile("Truc.txt");
  // for (int i = 0; i < p * (h+1); i++)
  //   OutputFile << V[i] << " ";
  // OutputFile.close();
  ClassesHeap H(V, int(p), int(h), int(NbCLasses));
  //  Rcpp::Rcout<<H;
  for(int index=0; index< 3*(p-1); index++) merge[index]=H.Output[index];
  return merge;
}

