
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector wcss_single(const arma::SpMat<double> & C, const NumericVector & cluster){

	NumericVector result(max(cluster));

	int baseline = 0;
	for(int i=1; i<cluster.length(); i++){
		if( cluster(i) != cluster(i-1) | i+1 == cluster.length() ){
			arma::SpMat<double> C_sub = C.submat(baseline, baseline, i, i);

			arma::sp_mat::iterator start = C_sub.begin();
			arma::sp_mat::iterator end = C_sub.end();

			double total = 0;
			for(arma::sp_mat::iterator it = start; it != end; ++it){
				total += (*it);
			}
			result[cluster(i-1)-1] = total;
			baseline = i;
		}
	}

	return result;
}
    

// [[Rcpp::export]]
NumericVector WCSS(const arma::SpMat<double> & C, const NumericMatrix & clusterMat){

	NumericVector result(clusterMat.ncol());

	#pragma omp parallel for if(parallelism_enabled)  
	for(int j=0; j<clusterMat.ncol(); j++){

		NumericVector res = wcss_single(C, clusterMat(_,j) );
		result[j] = sum(res);
	}

	return result;
}