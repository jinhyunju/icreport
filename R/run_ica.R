#' Automated ICA analysis with Report
#'
#' Performing ICA on a dataset and create a list object with results.
#'
#' @param phenotype.mx Phenotype matrix with diemnsions g x N   \code{phenotype.mx}
#' @param info.df Dataframe that holds sample covariates (ex. population, gender, age, etc...) \code{info.df}
#' @param check.covars Column names of info.df which hold the covariates
#' that should be used for association testing with IC coefficients \code{check.covars}
#' @param k.est Number of components to be estimated or method to estimate it. \code{k.est}
#' @param scale.pheno Logical value specifying the scaling of row of the phenotype.mx. \code{scale.pheno}
#' @param n.runs Number of runs for estimating k. Default value is set to 5. \code{n.runs}
#' @param max.iter Maximum iterations for estimating k for each run. Default value is set to 10. \code{max.iter}
#' @param n.cores Number of cores to be used for estimating k. Default is set to 1. \code{n.cores}
#' @param cor.threshold Threshold for significant correlation calling. Default is set to 0.05 \code{cor.threshold}
#' @return List with the following entries.
#' @keywords keywords
#'
#' @import mclust
#' @export
#'
#' @examples
#' R code here showing how your function works
run_ica <- function(phenotype.mx = NULL, info.df = NULL, check.covars = NULL,
                    k.est = NULL, scale.pheno = FALSE,
                    n.runs = 5, max.iter = 10, n.cores = 1, cor.threshold = 0.05, ...){

    if(is.null(phenotype.mx)){
        cat("Error: Phenotype matrix is missing \n")
        break;
    }
    # removing 0 variance genes and scaling and centering the phenotype matirx
    phenotype.mx <- pre_process_data(phenotype.mx, scale.pheno = scale.pheno)

    # Estimating the number of components
    k.est.plot.df <- NULL
    if(is.null(k.est)){
        cat("Method for estimating k = Default (limiting dimension of phenotype matrix) \n")

        k.est <- min(dim(phenotype.mx))

        cat("Number of Components = ",k.est,"\n")

    } else if (k.est == "repratio"){
        cat("Method for estimating k = Replicatio ratio \n")

        k.est.result <- k_est_multirun(phenotype.mx, n.runs = n.runs,
                                     max.iter = max.iter,n.core = n.cores)
        k.est <- k.est.result[["k.est"]]
        k.est.plot.df <- k.est.result[["plot.df"]]
        cat("Number of Components = ",k.est,"\n")
    }

    cat("Running ICA with ",k.est, " Independent Componenets \n")


    ica.result <- fastICA::fastICA(phenotype.mx, k.est,
                           alg.typ = "parallel",method = "C",
                           fun = "logcosh" ,                            # function that should be used to estimate ICs, default is logcosh
                           alpha = 1,row.norm = FALSE,                  # row.norm is set to false since the phenotype.mx is scaled separately
                           maxit=500,tol = 0.0001, verbose = FALSE)


#    ica.result <- fastICA_gene_expr(phenotype.mx, k.est)

    rownames(ica.result$S) <- rownames(phenotype.mx)                   # Setting appropriate names for signals and mixing matrix
    colnames(ica.result$A) <- colnames(phenotype.mx)
    rownames(ica.result$A) <- paste("IC",c(1:dim(ica.result$A)[1]),sep="")

    # Attaching the sample info dataframe to the ica list
    ica.result$info.df <- info.df

    if(!is.null(k.est.plot.df)){
        ica.result$k.plot.df <- k.est.plot.df
        rm(k.est.plot.df)
    }

    cat("Estimating Number of Peaks in each IC \n")
    ica.result$peak.results <- apply(ica.result$S, 2, peak_detection)
    # peaks are defined as gene contributions that are larger than 2 standard deviations
    ica.result$peak.mx <- apply(ica.result$S, 2, function(x) 1*(abs(x) > 2*sd(x)))

    cat("Calculating Variance Explained by each IC \n")
    # get the total variance by the sums of squares of the scaled phenotype.mx
    total.var <- sum(phenotype.mx^2)

    # applying IC component-wise variance calculations
    var.IC <- sapply(1:dim(ica.result$A)[1],
                     function (x) IC_variance_calc(ica.result$S[,x], ica.result$A[x,]))

    # % variance explained by each IC
    percent.var <- (var.IC / total.var) * 100

    cat("Sanity Check : Total % of variance explained by",k.est,"ICs = ", sum(percent.var), "\n")

    cat("Creating index based on Variance explained \n")
    ica.result$order <- order(percent.var,decreasing = T) # ordering the ICs based on the amount of variance they explain



    # Checking correlation between IC coefficients and measured covariates
    if(!is.null(check.covars)){
        # Anova analysis for covariates vs ICA weights (A matrix)
        ica.result$cov.pval.mx <- component_association_test(ica.result$A,info.df,check.covars)
        corr.idx <- which(ica.result$cov.pval.mx < cor.threshold, arr.ind = T)
    } else{
        corr.idx <- NULL
    }

    if(length(corr.idx) != 0 ){
        ica.result$cov.corr.idx <- data.frame("Signal.idx" = corr.idx[,1],                    # which IC is correlated with
                                              "Covariate.idx" = corr.idx[,2],                 # which covariate
                                              "Covariate.Name" = check.covars[corr.idx[,2]])  # with the name of
        rm(corr.idx)
    } else {
        ica.result$cov.corr.idx <- NULL     # in case there are no apparent covariates
    }



    if(is.null(ica.result$cov.corr.idx)){
        sig <- rep(0,k.est)
        correlated.ic <- NULL
    } else{
        sig <- rep(0,k.est)
        correlated.ic <- unique(ica.result$cov.corr.idx$Signal.idx)
        sig[correlated.ic] <- 1
    }

    mclust.result <- suppressMessages(apply(ica.result$A, 1, function(x) Mclust(x)))

    ica.result$ica.stat.df <- data.frame("N.peaks"=sapply(ica.result$peak.results, function(x) x$N), # Number of peaks for each IC
                                         "n.clust"= sapply(mclust.result, function(x) x$G),      # Number of predicted clusters
                                         "percent.var" = percent.var,                                     # Percent variance explained
                                         "corr.ic" = factor(sig), "idx" = c(1:k.est))             # if correlated with covariate = 1 , 0 otherwise


    # which IC has more than 1 predicted clusters?
    multi.clust <- which(ica.result$ica.stat.df$n.clust > 1)

    # get union between multi clusters and correlated ICs
    hf.vec <- sort(union(multi.clust,correlated.ic))

    hf.vec.names <-paste("IC",hf.vec,sep="")
    # name the ICs as IC#

    cat("ICs marked as confounding factors = \n")
    print(hf.vec.names)

    # Creating a matrix indicating which genes are influenced by a given IC
    ica.result$ica.confeti.mx <- matrix(0,nrow = dim(ica.result$S)[1], ncol = length(hf.vec))
    # rownames = gene names
    rownames(ica.result$ica.confeti.mx) <- rownames(ica.result$S)
    # column names = correlated ICs + multi cluster ICs
    colnames(ica.result$ica.confeti.mx) <- hf.vec.names

    cat("Creating CONFETI matrix for regression \n")
    for ( i in 1:length(hf.vec)){
      k <- hf.vec[i]
      ic.name <- hf.vec.names[i]
      peak.temp <- names(ica.result$peak.results[[k]]$peaks)
      ica.result$ica.confeti.mx[peak.temp,ic.name] <- 1
    }

    return(ica.result)
}
