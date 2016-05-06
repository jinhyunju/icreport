#' Generating ICA or PCA summary reports.
#'
#' Generating a HTML report from a ICA or PCA list object.
#'
#' @param ica_result ICA or PCA result list created by either \code{gene_expr_ica()} or \code{gene_expr_pca()}.
#' @param n.comps Number of principal components to plot.
#'        Default is set to plot every PC. Only used for PCA plotting not for ICA plotting.
#' @param prefix Output filename prefix. The output file will be named
#'        "prefix_ICA_summary.html".
#' @param output.path Directory path for generating the output HTML file.
#'        default is set to current working directory.
#' @param file.ext File extension to be used for saved plots.
#'        Default is set to png for html reports.
#'        Note that if you use pdf plots for html files they will not show up in your report.
#' @return output HTML report.
#' @keywords keywords
#'
#' @import ggplot2
#' @import knitr
#' @import rmarkdown
#'
#' @export
ica_report <- function(ica_result = NULL,
                       h5_file = NULL,
                       n.comps = NULL,
                       prefix = NULL,
                       output.path = NULL, file.ext = "png"){

    if(!exists("ica_result")){
        stop("Please specify the input to generate a report. \n")
    }

    if(!exists("prefix")){
       stop("Please specify the prefix of the output file \n")
    }

    if(is.null(output.path)){
        message("Output path is not specified, using current working directory \n")
        output.path <- getwd()

    }
    A_mx <- ica_result$A

    options(warn=-1)
    mclust.result <- apply(A_mx, 1, function(x) mclust::Mclust(x))
    # Turn warnings back on to prevent disasters
    options(warn=0)

    ica_summary_df <- data.frame("N.peaks"=sapply(ica_result$peaks, function(x) length(x)), # Number of peaks for each IC
                                 "percent.var" = ica_result$percent_var)                   # Percent variance explained
    ica_summary_df$idx <- 1:nrow(ica_summary_df)

    threshold <- 0.05

    geno_mx <- load_h5_data(h5_file, "genotypes")
    # change condition to something more general
    if(attr(ica_result, "geno_cor") == "yes"){
      sig <- rep(0,nrow(A_mx))
      sig_geno <- rep(0,nrow(A_mx))

      ica_summary_df$geno_cor <- ica_result$geno_cor_count
      ica_summary_df$geno_cor_max <- apply(ica_result$geno_pvals, 2, function(x) which.min(x))
    }

    # change condition to something more general
    if(attr(ica_result, "covar_cor") == "yes" ){
      covar_sig_idx <- which(ica_result$covar_pvals < (threshold / (length(ica_result$covar_pvals))), arr.ind = T)

      sig <- rep(0,nrow(A_mx))
      sig_covars <- rep(0,nrow(A_mx))

      if(length(covar_sig_idx) != 0 ){
        message(nrow(covar_sig_idx ), " associations detected")
        correlated.ic <- unique(covar_sig_idx [,1])

        for( c in 1:length(correlated.ic)){
          ic.index <- correlated.ic[c]
          sig_covars[ic.index] <- which.min(ica_result$covar_pvals[ic.index,])
          sig[ic.index] <- 1
        }
      } else {
        message("No ICs were correlated with covariates")

      }
      ica_summary_df$covar_cor <- sig
      ica_summary_df$covar_cor_idx <- sig_covars


    }

    if((attr(ica_result, "geno_cor") == "yes") & (attr(ica_result, "covar_cor") == "yes")){
      ica_summary_df <- ica_summary_df[with(ica_summary_df, order(geno_cor, covar_cor ,percent.var, decreasing = TRUE)),]
    } else if ((attr(ica_result, "geno_cor") == "yes") & (attr(ica_result, "covar_cor") == "no")) {
      ica_summary_df <- ica_summary_df[with(ica_summary_df, order(geno_cor,percent.var, decreasing = TRUE)),]
    } else if ((attr(ica_result, "geno_cor") == "no") & (attr(ica_result, "covar_cor") == "yes")) {
      ica_summary_df <- ica_summary_df[with(ica_summary_df, order(covar_cor ,percent.var, decreasing = TRUE)),]
    } else {
      ica_summary_df <- ica_summary_df[with(ica_summary_df, order(percent.var, decreasing = TRUE)),]
    }


    # create axis tick positions for chromosomes if gene info is supplied
    # the same axis tick positions are used for every plot so it only needs to be calculated once
    x.axis <- chr_axis_creator(ica_result$gene_info)

    markdown.file <- system.file("templates/ICA_Visualization.Rmd", package="icreport")

    outFile = paste(output.path,"/",prefix,"_ICA_summary.html",sep="")

 #   suppressMessages(rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document"))

    knitr::knit2html(markdown.file,output = outFile)



}
