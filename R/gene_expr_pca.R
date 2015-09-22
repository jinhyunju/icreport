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
gene_expr_pca <- function(phenotype.mx = NULL, info.df = NULL, check.covars = NULL,
                          scale.pheno = FALSE, cor.threshold = 0.05){

    if(is.null("phenotype.mx")){
      message("Please specify phenotype matrix \n")
      break;
    }

    pheno.nrow <- nrow(phenotype.mx)
    pheno.ncol <- ncol(phenotype.mx)
    message("phenotype.mx dimensions = ", pheno.nrow , " x ", pheno.ncol, "\n")

    if (pheno.nrow < pheno.ncol){
      message("[Caution] Number of samples exceeding number of measured features,
              please check rows and columns of <phenotype.mx> \n")
      message("- If you are from the future and have more samples than measured features,
              disregard the above message and please proceed. \n")
    }

    if(!is.null(info.df)){

      message("Checking dimensions and column/row names \n")

      if(nrow(info.df) == ncol(phenotype.mx)){
        message("[Good] Number of samples in <info.df> and <phenotype.mx> match \n")
      } else {
        stop("Error: Sample numbers in <info.df> and <phenotype.mx> don't match. Stopping script")
      }

      matching.names <- sum(rownames(covars) %in% colnames(phenotype.mx))

      if(matching.names == ncol(phenotype.mx)){
        message("[Good] All samples in <phenotype.mx> are accounted for in <info.df> \n")
      } else {
        stop("Missing sample Information: Check rownames of <info.df> and column names of <phenotype.mx>")
      }
    }

    if(!is.null(info.df) & is.null(check.covars)){
      message("- Sample info detected but missing input for <check.covars> \n")
      message("- Using column names of <info.df> as <check.covars> \n")
      check.covars <- colnames(info.df)
    }

    message("Pre Processing Data \n")
    # removing NA values, centering, scaling (optional)
    phenotype.mx <- pre_process_data(phenotype.mx,
                                     scale.pheno = scale.pheno)

    message("Running PCA \n")
    pca.result <- prcomp(t(phenotype.mx))

    # calculating percent variance each component explains
    pca.result$var.percent <- (pca.result$sdev^2 / sum(pca.result$sdev^2)) * 100


    pca.result$peaks <- apply(pca.result$rotation, 2, peak_detection)

    if(!is.null(check.covars) & !is.null(info.df)){
      message("Checking association between covariates and components\n")
      # Anova analysis for covariates vs ICA weights (A matrix)
      pca.result$cov.pval.mx <- component_association_test(t(pca.result$x),info.df,check.covars)

      # get bonferroni corrected threshold
      cor.threshold <- cor.threshold / (length(pca.result$cov.pval.mx))

      # row and column positions for significant p-values
      corr.idx <- which(pca.result$cov.pval.mx < cor.threshold, arr.ind = T)
    } else{
      message("No <info.df> supplied, association test skipped.\n")
      corr.idx <- NULL
    }

    if(length(corr.idx) != 0 ){

      # get PCs that were correlated to a covariate
      correlated.pc <- unique(corr.idx[,1])

      # Initiate dataframe to save results
      covariate.corr.df <- data.frame(matrix(nrow = length(correlated.pc),ncol = 3))
      colnames(covariate.corr.df) <- c("PC","Covariate.idx","Covariate.Name")

      for( c in 1:length(correlated.pc)){
        # get the PC number
        pc.index <- correlated.pc[c]

        # get the most significantly correlated covariate
        covar.index <- which.min(pca.result$cov.pval.mx[pc.index,])
        covariate.corr.df[c,"PC"] <- pc.index
        covariate.corr.df[c,"Covariate.idx"] <- covar.index
        covariate.corr.df[c,"Covariate.Name"] <- names(covar.index)
      }

      pca.result$cov.corr.idx <- covariate.corr.df

      # Add variance information
      pca.result$cov.corr.idx$var <- pca.result$var.percent[pca.result$cov.corr.idx$PC]

      # remove temprary objects just for clarity
      rm(covariate.corr.df, covar.index, pc.index, c)
    } else {
      pca.result$cov.corr.idx <- NULL     # in case there are no associated covariates
    }

    # save sample information in output object
    pca.result$info.df <- info.df

    # label result as pca for report2me() function
    attr(pca.result, 'method') <- "pca"
    return(pca.result)
}
