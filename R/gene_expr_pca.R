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


    pca.result$var.percent <- pca.result$sdev / sum(pca.result$sdev) * 100

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
      correlated.pc <- unique(corr.idx[,1])

      covariate.corr.df <- data.frame(matrix(nrow = length(correlated.pc),ncol = 3))
      colnames(covariate.corr.df) <- c("PC","Covariate.idx","Covariate.Name")
      for( c in 1:length(correlated.pc)){
        pc.index <- correlated.pc[c]
        covar.index <- which.min(pca.result$cov.pval.mx[pc.index,])
        covariate.corr.df[c,"PC"] <- pc.index
        covariate.corr.df[c,"Covariate.idx"] <- covar.index
        covariate.corr.df[c,"Covariate.Name"] <- names(covar.index)
      }

      pca.result$cov.corr.idx <- covariate.corr.df
      pca.result$cov.corr.idx$var <- pca.result$var.percent[pca.result$cov.corr.idx$PC]

      rm(covariate.corr.df, covar.index, pc.index, c)
    } else {
      pca.result$cov.corr.idx <- NULL     # in case there are no associated covariates
    }


    pca.result$info.df <- info.df

    return(pca.result)
}
