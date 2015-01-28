#' Custom PCA function for analyzing gene expression data.
#'
#' Performing PCA on a dataset and create a list object with results.
#'
#' @param phenotype.mx Phenotype matrix with diemnsions g x N
#' @param info.df Dataframe that holds sample covariates (ex. population, gender, age, etc...)
#' @param check.covars Column names of info.df which hold the covariates
#' that should be used for association testing with IC coefficients.
#' @param scale.pheno Logical value specifying the scaling of row of the phenotype.mx. Default is set to FALSE.
#' @param cor.threshold Threshold for significant correlation calling. Default is set to 0.05.
#' @return List with the following entries.
#' @keywords keywords
#'
#' @import mclust
#' @export
#'
#' @examples
#' R code here showing how your function works
gene_expr_pca <- function(phenotype.mx = NULL, info.df = NULL, check.covars = NULL,
                          scale.pheno = FALSE, cor.threshold = 0.05, ...){

    if(is.null("phenotype.mx")){
      cat("Please specify phenotype matrix \n")
      break;
    }

    cat("Running PCA\n")
    pca.result <- prcomp(t(phenotype.mx))

    if(!is.null(check.covars) & !is.null(info.df)){
      cat("Checking association between covariates and components\n")
      # Anova analysis for covariates vs ICA weights (A matrix)
      pca.result$cov.pval.mx <- component_association_test(t(pca.result$x),info.df,check.covars)
      cor.threshold <- cor.threshold / (dim(pca.result$cov.pval.mx)[1] * dim(pca.result$cov.pval.mx)[2])
      corr.idx <- which(pca.result$cov.pval.mx < cor.threshold, arr.ind = T)
    } else{
      cat("No info.df supplied, association test skipped.\n")
      corr.idx <- NULL
    }

    if(length(corr.idx) != 0 ){
      pca.result$cov.corr.idx <- data.frame("Signal.idx" = corr.idx[,1],                    # which IC is correlated with
                                            "Covariate.idx" = corr.idx[,2],                 # which covariate
                                            "Covariate.Name" = check.covars[corr.idx[,2]])  # with the name of
      rm(corr.idx)
    } else {
      pca.result$cov.corr.idx <- NULL     # in case there are no associated covariates
    }
    pca.result$info.df <- info.df

    return(pca.result)
}
