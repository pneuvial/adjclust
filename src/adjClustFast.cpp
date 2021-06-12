// #define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#if defined(_OPENMP)
#include <omp.h>
extern const bool parallelism_enabled = true;
#else
extern const bool parallelism_enabled = false;
#endif

// write RcppExports.R with: 
//       Rcpp::compileAttributes(".")



// Exports from mainFunctions.c

// There are written to RcppExport.cpp
// 	Once they are there, comment this out

// // [[Rcpp::export]]
// extern "C" SEXP percDown(SEXP Rpositions, SEXP Rdistances, SEXP Rl, SEXP Rpos);

// // [[Rcpp::export]]
// extern "C" SEXP cWardHeaps(SEXP RrcCumRight, SEXP RrcCumLeft, SEXP Rh, SEXP Rp, SEXP RchainedL, SEXP Rpositions, SEXP Rdistances, SEXP RlHeap, SEXP Rmerge, SEXP Rgains, SEXP RtraceW);




// Replace adjclust::matR() and adjclust::matL() with RcppArmadillo code
// Use separate version for matrix and sparseMatrix

// [[Rcpp::export]]
arma::SpMat<double> matL_sparse(const arma::SpMat<double> & Csq, const int & h) {

	int p = Csq.n_rows;
	arma::SpMat<double> out(p, h+1);

	// #pragma omp parallel for if(parallelism_enabled)              
	for( int i=0; i<p; i++){

		int k = 0;
		double value;

		for( int j=i; j<std::min(i+h+1,p); j++){

			value = Csq(i,j);

			if(k == 0){
		 		out(i,k) = value;
			}else{
				if( value!= 0.0){
					out(i,k) = 2.0*value;
				}
			}
			k++;
		}
	}

	return out;	           
}


// [[Rcpp::export]]
arma::mat matL_full(const arma::mat & Csq, const int & h) {

	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

	// #pragma omp parallel for if(parallelism_enabled)                  
	for( int i=0; i<p; i++){
		
		int k = 0;
		double value;

		for( int j=i; j<std::min(i+h+1,p); j++){

			value = Csq(i,j);

			if(k == 0){
		 		out(i,k) = value;
			}else{
				out(i,k) = 2.0*value;
			}
			k++;
		}
	}

	return out;	           
}




// [[Rcpp::export]]
arma::SpMat<double> matR_sparse(const arma::SpMat<double> & Csq, const int & h) {
	               
	int p = Csq.n_rows;
	arma::SpMat<double> out(p, h+1);

	// #pragma omp parallel for if(parallelism_enabled)    
	for( int i=0; i<p; i++){
		
		int k = 0;
		double value;

		for( int j=i; j>=std::max(i-h, (int) 0); j--){
		   value = Csq(i,j);
		  if(k == 0){
		    out(p-i-1,k) = value;
		  }else{
		    if( value != 0.0){
		        out(p-i-1,k) = 2.0*value;
		    }
		  }
		  k++;
		}
	}

	return out;
}


// [[Rcpp::export]]
arma::mat matR_full(const arma::mat & Csq, const int & h) {
	               
	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

	// #pragma omp parallel for if(parallelism_enabled)    
	for( int i=0; i<p; i++){
		
		int k = 0;
		double value;
		
		for( int j=i; j>=std::max(i-h, (int) 0); j--){

		  value = Csq(i,j);

		  if(k == 0){
		    out(p-i-1,k) = value;
		  }else{
		  	out(p-i-1,k) = 2.0*value;
		  }
		  k++;
		}
	}

	return out;
}
