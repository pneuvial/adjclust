// #define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#if defined(_OPENMP)
  #include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

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
arma::mat matL_sparse_rowCumsums(const arma::SpMat<double> & Csq, const int & h, int nthreads) {

	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(nthreads)
  #endif
	
	for( int i=0; i<p; i++){
		
		int k = 0;
		double value;

		for( int j=i; j<std::min(i+h+1,p); j++){

			value = Csq(i,j);

			if(k == 0){
		 		out(i,k) = value;
			}else{
				out(i,k) = out(i,k-1) + 2.0*value;
			}
			k++;
		}

		// finish cumsum
		while(k < h+1){			
		  	out(i,k) = out(i,k-1);
		  	k++;
		}
	}

	return out;	           
}




// [[Rcpp::export]]
arma::mat matL_full_rowCumsums(const arma::mat & Csq, const int & h, int nthreads) {

	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(nthreads)
  #endif

	for( int i=0; i<p; i++){
		
		int k = 0;
		double value;

		for( int j=i; j<std::min(i+h+1,p); j++){

			value = Csq(i,j);

			if(k == 0){
		 		out(i,k) = value;
			}else{
				out(i,k) = out(i,k-1) + 2.0*value;
			}
			k++;
		}

		// finish cumsum
		while(k < h+1){			
		  	out(i,k) = out(i,k-1);
		  	k++;
		}
	}

	return out;	           
}

// [[Rcpp::export]]
arma::SpMat<double> matR_sparse(const arma::SpMat<double> & Csq, const int & h) {
	               
	int p = Csq.n_rows;
	arma::SpMat<double> out(p, h+1);

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

// [[Rcpp::export]]
arma::mat matR_sparse_rowCumsums(const arma::SpMat<double> & Csq, const int & h, int nthreads) {
	               
	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(nthreads)
  #endif

	for( int i=0; i<p; i++){
		
		int k = 0;
		double value;
		
		for( int j=i; j>=std::max(i-h, (int) 0); j--){

		  value = Csq(i,j);

		  if(k == 0){
		    out(p-i-1,k) = value;
		  }else{
		  	// set value while computing cumulative row sum
		  	out(p-i-1,k) = out(p-i-1,k-1) + 2.0*value;
		  }
		  k++;
		}

		// finish cumsum
		while(k < h+1){			
		  	out(p-i-1,k) = out(p-i-1,k-1);
		  	k++;
		}
	}

	return out;
}

// [[Rcpp::export]]
arma::mat matR_full_rowCumsums(const arma::mat & Csq, const int & h, int nthreads) {
	               
	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(nthreads)
  #endif

	for( int i=0; i<p; i++){
		
		int k = 0;
		double value;
		
		for( int j=i; j>=std::max(i-h, (int) 0); j--){

		  value = Csq(i,j);

		  if(k == 0){
		    out(p-i-1,k) = value;
		  }else{
		  	// set value while computing cumulative row sum
		  	out(p-i-1,k) = out(p-i-1,k-1) + 2.0*value;
		  }
		  k++;
		}

		// finish cumsum
		while(k < h+1){			
		  	out(p-i-1,k) = out(p-i-1,k-1);
		  	k++;
		}
	}

	return out;
}
