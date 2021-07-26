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
arma::mat matL_sparse_rowCumsums(const arma::SpMat<double> & Csq, const int & h) {

	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

	#pragma omp parallel for if(parallelism_enabled)                  
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
arma::mat matL_full_rowCumsums(const arma::mat & Csq, const int & h) {

	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

	#pragma omp parallel for if(parallelism_enabled)                  
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



// [[Rcpp::export]]
arma::mat matR_sparse_rowCumsums(const arma::SpMat<double> & Csq, const int & h) {
	               
	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

	#pragma omp parallel for if(parallelism_enabled)    
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
arma::mat matR_full_rowCumsums(const arma::mat & Csq, const int & h) {
	               
	int p = Csq.n_rows;
	arma::mat out(p, h+1, arma::fill::zeros);

	#pragma omp parallel for if(parallelism_enabled)    
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




using namespace Rcpp;

// [[Rcpp::export]]
NumericVector wcss_single(const arma::SpMat<double> & C, const NumericVector & cluster){

	NumericVector result(max(cluster));
	arma::SpMat<double> C_sub;

	int baseline = 0;
	double n_features;

	// for each feature 
	for(int i=0; i<cluster.length(); i++){

		// Split conditions 
		// if the first condition is TRUE, then cluster(i+1) causes overflow
		bool flag=FALSE;
		if( i+1 == cluster.length() ){
			flag = TRUE;
		}else if( cluster(i) != cluster(i+1) ){
			flag = TRUE;
		}

		// if current and previous features belong to different features
		if( flag ){

			double total = 0;
			if( baseline == i){
				total = C(i,i);
			}else{
				// get submatrix corrsponding to this cluster
				C_sub = C.submat(baseline, baseline, i, i);

				// Evalute sum of all matrix elements
				for(arma::sp_mat::iterator it = C_sub.begin(); it != C_sub.end(); ++it){
					total += (*it);
				}
			}
			// total sum of squares divided by number of features in the group
			n_features = (i - baseline + 1);
			result[cluster(i)-1] = (total - n_features) / n_features;
			baseline = i+1;
		}
	}

	return result;
}
    

// [[Rcpp::export]]
NumericVector WCSS(const arma::SpMat<double> & C, const NumericMatrix & clusterMat){

	std::vector<double> result(clusterMat.ncol());

	#pragma omp parallel for if(parallelism_enabled)
	for(int j=0; j<clusterMat.ncol(); j++){

		NumericVector v = clusterMat(_,j);

		// for each number of clusters, compute within-cluster sum of squares
		double Total = sum( wcss_single(C, v) );

		result[j] = Total;
	}

	return wrap(result);
}


