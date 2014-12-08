#' Estimating the number of components (k)
#'
#' Estimating the number of components to be estimated
#' through the replication ratio between two independent runs of ICA
#'
#' @param phenotype.mx The matrix that holds the gene expression profile   \code{phenotype.mx}
#' @param max.iter Number of maximum iterations, default is set to 10 \code{max.iter}
#' @param run index \code{run.index}
#' @return output A matrix holding the p-values for each indepenedent component and covariate pair.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
k_est_multirun <- function(phenotype.mx, n.runs = 5, max.iter = 10,n.cores = 1){

    if( n.cores == 1){
        temp.result <- lapply(1:n.runs, function(x) k_est_singlerun(phenotype.mx,
                                                               max.iter = max.iter,run = x))

    } else if (n.cores > 1) {
        # parallel version of the code above
        temp.result <- parallel::mclapply(1:n.runs,
                                          function(x) k_est_singlerun(phenotype.mx, max.iter = max.iter,run = x),
                                          mc.cores = n.cores)
    }

    # get all k estimates from multiple runs
    k.vec <- sapply(temp.result, function(x) x$k.est)

    # get the summary data for each estimation run
    summary.df <- Reduce(function(...) merge(..., all=T),
                       lapply(temp.result, function(x) x$summary.vec ))

    k.est <- max(k.vec)
    cat("Final estimate of k = ",k.est,"\n")

    # saving the result in a dataframe
    plot.df <- data.frame(summary.df)
    plot.df$run <- as.factor(plot.df$run)


    return(list("k.est" = k.est, "plot.df" = plot.df))
}
