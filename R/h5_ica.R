#' Custom ICA function for analyzing gene expression data.
#'
#' Performing ICA on a dataset and create a list object with results.
#' @param h5_file Name of HDF5 file that has the phenotype values saved.
#' @param input_pheno Phenotype matrix with diemnsions N x g
#' @param k.est Number of components to be estimated or method to estimate it.
#' @param scale.pheno Logical value specifying the scaling of row of the phenotype.mx.
#' @return List with the following entries.
#' @keywords keywords
#'
#' @export
h5_ica <- function(h5_file = NULL,
                   input_pheno = NULL,
                   k.est = NULL,
                   scale.pheno = FALSE,
                   var.cutoff = 99){

    if(is.null(h5_file) & is.null(input_pheno)){
        stop("Error: Input HDF5 file missing")
    } else if (is.null(h5_file) & !is.null(input_pheno)){
        message("Using in-memory input object as phenotype matrix \n")
        phenotype.mx <- input_pheno
    } else if (!is.null(h5_file)){
        message("Loading phenotype values from hdf5 file = ", h5_file,"\n")
        phenotype.mx <- load_h5_data(h5_file, "phenotypes")

    }

    if(is.null(colnames(phenotype.mx))){
      message("<phenotype.mx> is missing column names, set to default. \n")
      colnames(phenotype.mx) <- paste("sample",c(1:ncol(phenotype.mx)), sep = "_")
    }

    if(is.null(rownames(phenotype.mx))){
      message("<Phenotype.mx> is missing row names, set to default. \n")
      rownames(phenotype.mx) <- paste("feature",c(1:nrow(phenotype.mx)), sep = "_")
    }

    pheno.nrow <- nrow(phenotype.mx)
    pheno.ncol <- ncol(phenotype.mx)

    message("Original dimensions of <phenotype.mx> = ", pheno.nrow , " x ", pheno.ncol, "\n")

    if (pheno.nrow > pheno.ncol){
      message("[Caution] Number of samples exceeding number of measured features,
              please check rows and columns of <phenotype.mx> \n")
      message("- If you are from the future and have more samples than measured features,
              disregard the above message. \n")
    }

    message("------ Pre-processing Data ------- \n")
    # removing 0 variance genes and scaling and centering the phenotype matirx
    phenotype.mx <- pre_process_data(t(phenotype.mx), scale.pheno = scale.pheno)

    if(is.null(k.est)){

      message("Missing input for <k.est>, using the number of principal components explaining ", var.cutoff, "% of total variance \n")
      pca.pheno <- prcomp(t(phenotype.mx))
      percent <- (cumsum(pca.pheno$sdev^2) /sum(pca.pheno$sdev^2)) * 100
      k.est <- which(percent > var.cutoff)[1]
      message(k.est," components needed to explain more than ",var.cutoff,"% of the variance \n")

      if(k.est == 1){
        stop("1 component explains more than",var.cutoff,"% of the variance,
             check your data or set <k.est> to a number bigger than 1 \n")
      }

    }
    message("Running ICA with ", k.est," components to be estimated")
    message("* This may take some time depending on the size of your dataset \n")

    ica.result <- fastICA_gene_expr(phenotype.mx, k.est,
                                    fun = "logcosh",                            # function that should be used to estimate ICs, default is logcosh
                                    alpha = 1, scale.pheno = FALSE,                  # row.norm is set to false since the phenotype.mx is scaled separately
                                    maxit=500, tol = 0.0001, verbose = FALSE)

    rownames(ica.result$S) <- rownames(phenotype.mx)                   # Setting appropriate names for signals and mixing matrix
    colnames(ica.result$S) <- paste("IC",c(1:dim(ica.result$S)[2]),sep="")
    colnames(ica.result$A) <- colnames(phenotype.mx)
    rownames(ica.result$A) <- paste("IC",c(1:dim(ica.result$A)[1]),sep="")
    message("------ Post-ICA Processing ------- \n")
    message("- Labeling peak genes in each IC \n")
    ica.result$peaks <- apply(ica.result$S, 2, peak_detection)
    # peaks are defined as gene contributions that are larger than 2 standard deviations
    ica.result$peak.mx <- apply(ica.result$S, 2, function(x) 1*(abs(x) > 2*sd(x)))

    message("- Calculating variance explained by each IC \n")
    # get the total variance by the sums of squares of the scaled phenotype.mx
    total.var <- sum(phenotype.mx^2)

    # applying IC component-wise variance calculations
    var.IC <- sapply(1:dim(ica.result$A)[1],
                     function (x) IC_variance_calc(ica.result$S[,x], ica.result$A[x,]))

    percent.var <- (var.IC / total.var) * 100
    ica.result$percent_var <- percent.var
    message("- Sanity Check : Total % of variance explained by ",k.est," ICs = ", sum(percent.var), "\n")

    class(ica.result) <- "ICAobject"
    attr(ica.result, 'method') <- "ica"
    attr(ica.result, 'covar_cor') <- "no"
    attr(ica.result, 'geno_cor') <- "no"
    return(ica.result)
}

#' Custom ICA function for analyzing gene expression data.
#'
#' Performing ICA on a dataset and create a list object with results.
#' @param h5_file Name of HDF5 file that has the phenotype values saved.
#' @param input_pheno Phenotype matrix with diemnsions N x g
#' @param k.est Number of components to be estimated or method to estimate it.
#' @param scale.pheno Logical value specifying the scaling of row of the phenotype.mx.
#' @return List with the following entries.
#' @keywords keywords
#'
#' @import mclust
#' @export
ica_covar_check <- function(ica_list = NULL,
                            h5_file = NULL,
                            covars = NULL, cor.threshold = 0.05){

    if(is.null(covars) & !is.null(h5_file)){
      message("Loading covariates from HDF5 file \n")
      covars <- load_h5_data(h5_file, "covars")
    }

    A_mx <- ica_list$A

    message("Checking dimensions and column/row names \n")

    if(nrow(covars) == ncol(A_mx)){
      message("[Good] Number of samples in <covars> and <ica_list> match \n")
    } else {
      stop("Error: Sample numbers in <covars> and <ica_list> don't match. Stopping script")
    }

    matching.names <- sum(rownames(covars) %in% colnames(A_mx))

    if(matching.names == nrow(covars)){
      message("[Good] All samples in <ica_list> are accounted for in <covars> \n")
    } else {
      stop("Missing sample Information: Check sample names in <covars> and <ica_list>")
    }

    message("- Checking associations between ICs and covariates \n")
    # Anova analysis for covariates vs ICA weights (A matrix)
    cov.pval.mx <- ic_covariate_association_test(A_mx,covars)
    ica_list$covar_pvals <- cov.pval.mx
    attr(ica_list, 'covar_cor') <- "yes"

    return(ica_list)
}


#' Testing the association between independent component coefficients and known covariates.
#'
#' The function takes in the A matrix from the ica.result object
#' with a dataframe of measured covariates to test the association
#' between them.
#'
#' @param input.A A matrix of an ica_result list
#' @param info.input A dataframe that holds the measured covariates for each sample
#' @return output A matrix holding the p-values for each indepenedent component and covariate pair.
#' @keywords keywords
#'
#' @export
#'
ic_covariate_association_test <- function(input.A, info.input){

  n.components <- dim(input.A)[1]
  covar.names <- colnames(info.input)
  n.covariates <- length(covar.names)

  pval.mx <- matrix(NA, nrow = n.components, ncol = n.covariates)
  colnames(pval.mx) <- covar.names

  covar.df <- data.frame(info.input[colnames(input.A),covar.names])

  for( i in 1: dim(input.A)[1]){
    for (j in 1:length(covar.names)){
      analysis.df <- data.frame("IC"=input.A[i,],"Covar" = covar.df[,j])
      anova.fit <- aov(IC ~ Covar, data=analysis.df)
      model.fstat <- summary.lm(anova.fit)$fstatistic
      p.val <- pf(model.fstat[1],model.fstat[2],model.fstat[3], lower.tail = F, log.p = F)
      pval.mx[i,j] <- p.val
    }
  }

  return(pval.mx)
}

#' @import lrgpr
#' @import formula.tools
#' @export
ica_genotype_association <- function(ica_list = NULL,
                                     h5_file = NULL,
                                     genotype.mx = NULL,
                                     n.cores = 1){
  if(is.null(ica_list)){
      stop("ICA input missing")
  }

  if(is.null(genotype.mx) & !is.null(h5_file)){
    message("Loading genotypes from HDF5 file \n")
    genotype.mx <- load_h5_data(h5_file, "genotypes")
  } else if (is.null(genotype.mx) & is.null(h5_file)){
    stop("Missing input for genotypes, please specify genotype object or h5 file")
  }

  ica.loadings <- t(ica_list$A)

  ic.vs.geno <- lrgpr::glmApply(ica.loadings ~ SNP,
                         features = genotype.mx,
                         nthreads = n.cores)$pValues

  colnames(ic.vs.geno) <- rownames(ica_list$A)
  #sig <- which(ic.vs.geno < (0.05/length(ic.vs.geno) ), arr.ind = TRUE)

  #genetic.factors <- colnames(ic.vs.geno)[unique(sig[,"col"])]
  #non.genetic <- colnames(ic.vs.geno)[which(!(colnames(ic.vs.geno) %in% genetic.factors))]
  ica_list$geno_pvals <- ic.vs.geno
  attr(ica_list, 'geno_cor') <- "yes"
  return(ica_list)
}

#' @export
get_gene_info <- function(ica_list, h5_file){
#  gene_id <- h5read(h5_file, "phenotypes/col_info/id")
#  gene_chr <- h5read(h5_file, "phenotypes/col_info/pheno_chr")
#  gene_start <- h5read(h5_file, "phenotypes/col_info/pheno_start")
  gene_info_df <- as.data.frame(h5read(h5_file, "phenotypes/col_info"), stringsAsFactors = FALSE)
  gene_info_df$idx <- 1:nrow(gene_info_df)
  ica_list$gene_info <- gene_info_df
  return(ica_list)
}

#' @export
as.character.formula <- function(x){
  Reduce( paste, deparse(x) )
}
