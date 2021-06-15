# Gabriel Hoffman
# June 14, 2021
#
# Compute within cluster sum of squares

#' Compute within-clstuer sum of squares
#' 
#' Compute within-clstuer sum of squares for many values of cluster number.
#'
#' @param hcl hierarchical clustering result, where \code{hcl$data} is a similarity matrix
#' @param k_array array storing the number of clusters where the GAP statistic is evaluated.  All missing integer within \code{range(k_array)} are estimated using a linear interpretation of the GAP statistic.
#' @param p.value p-value cutoff used to determine the number of clusters
#' @param mc.cores number of CPU cores to use.
#' 
#' @return
#' \itemize{
#'   \item df_approx - data.frame storing results for each value of k
#'   \item n_clusters - number of clusters using p-value cutoff
#' }
#'
#' @details
#' The Within-Cluster Sum of Squares is computed for many values of cluster number, k.    Since evaluating the WCSS for each k is expensive, the WCSS for missing values of k are filled in using a linear interpolation. WCSS is used as to evaluted the number of clusters in the data, and is part of the GAP score \insertCite{tibshirani2001estimating}{adjclust}.  Here we adapted the score to be used for similarity matrices \insertCite{dehman2015performance}{adjclust}.
#'
#' @references
#' \insertRef{dehman2015performance}{adjclust}
#'
#' \insertRef{tibshirani2001estimating}{adjclust}
#'
#' @importFrom MASS fitdistr
#' @importFrom stats approxfun cutree
#' @importFrom stats pgamma var
#' @importFrom Matrix sparse.model.matrix crossprod diag
#' @import Rdpack
#' @importFrom parallel mclapply
#' @export
wcss = function(hcl, k_array, p.value=0.001, mc.cores=1){

	# make sure k values are unique and sorted
	k_array = sort(unique(k_array))

	# compute multiple cuts at the same time
	# this increase memory usage, but sames substantial time
	clustersMat = cutree(hcl, k=k_array)
	
	# Evaluate within-cluster sum of squares
	W = WCSS(hcl$data, clustersMat)

	df = data.frame(k = k_array, W = W)

	if( !any(df$W > 0) ){
		stop("Cannot interpolate when all W values are 0.")
	}

	# since evaluating each value of k is expensice,
	# interpolate for values of k that are skipped 
	f = approxfun(df$k, df$W)
	kmax = max(df$k)
	df_approx = data.frame(k = 1:kmax, W = f(1:kmax))
	df_approx$d = c(0,-diff(df_approx$W) / diff(df_approx$k))
	df_approx$Interpolated = TRUE
	df_approx$Interpolated[df_approx$k %in% k_array] = FALSE

	# fit gamma distribution to differences
	d = df_approx$d[-1]
	d = d[d>0] # only keep values greater than zero
	momAlpha <- (mean(d)^2)/var(d)
	momBeta <- var(d)/mean(d)
	suppressWarnings({
	fit = fitdistr(d, "gamma", start=list(shape=momAlpha, rate=momBeta))$estimate 
	})

	# evalute p-values based on the distribution of d values
	# assume that most d values are drawn from the null
	df_approx$pValues = 0
	df_approx$pValues[-1] = pgamma(d, shape=fit[1], rate=fit[2], lower.tail=FALSE)

	# Compute the last cluster with p-value less than the cutoff
	n_clusters = max(cumsum(df_approx$pValues < p.value))

	list(df_approx = df_approx, n_clusters = n_clusters)
}
