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
                                    alpha = 1,           # row.norm is set to false since the phenotype.mx is scaled separately
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
