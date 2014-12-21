#' Automated PCA report generator
#'
#' Generating a HTML report from a PCA object created with the function prcomp().
#'
#' @param pca.result PCA result list created by the function prcomp().   \code{pca.result}
#' @param info.df Sample information dataframe to be used for covariate correlation \code{info.df}
#' @param check.covars Names of covariates in info.df (column names)
#' that should be used for association tests \code{check.covars}
#' @param plot.pc Number of principal components to plot. Default is set to plot every PC. \code{plot.pc}
#' @param cor.threshold Threshold in calling an association between a covariate and a PC signiciant.
#'        Default value is set to 0.05. \code{cor.threshold}
#' @param prefix Output filename prefix. The output file will be named
#'        "prefix_ICA_summary.html". \code{prefix}
#' @param output.path Directory path for generating the output HTML file.
#'        default is set to current working directory. \code{output.path}
#'
#' @param file.ext Extension of plots saved in the process of compiling the HTML report.
#'        Default is set to png format. \code{file.ext}
#' @return output HTML report that contains component wise visualization and summary plots.
#' @keywords keywords
#'
#' @import ggplot2
#' @import mclust
#' @import knitr
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
pc_report <- function(phenotype.mx = NULL, info.df = NULL, check.covars = NULL, plot.pc = NULL,
                      cor.threshold = 0.05, prefix = NULL, output.path = NULL, file.ext = "png"){

    if(is.null("phenotype.mx")){
        cat("Please specify phenotype matrix \n")
        break;
    }

    if(is.null("prefix")){
        cat("Please specify the prefix of the output file \n")
        break;
    }

    if(is.null(output.path)){
        cat("Output path is not specified, using current working directory \n")
        output.path <- getwd()

    }

    pca.result <- prcomp(t(phenotype.mx))

    if(is.null(plot.pc)){
        plot.pc <- dim(pca.result$x)[2]
    }

    if(!is.null(check.covars) & !is.null(info.df)){
        # Anova analysis for covariates vs ICA weights (A matrix)
        pca.result$cov.pval.mx <- component_association_test(t(pca.result$x),info.df,check.covars)
        corr.idx <- which(pca.result$cov.pval.mx < cor.threshold, arr.ind = T)
    } else{
        corr.idx <- NULL
    }

    if(length(corr.idx) != 0 ){
        pca.result$cov.corr.idx <- data.frame("Signal.idx" = corr.idx[,1],                    # which IC is correlated with
                                              "Covariate.idx" = corr.idx[,2],                 # which covariate
                                              "Covariate.Name" = check.covars[corr.idx[,2]])  # with the name of
        rm(corr.idx)
    } else {
        pca.result$cov.corr.idx <- NULL     # in case there are no apparent covariates
    }


    data.set <- prefix
    markdown.file <- system.file("templates/PCA_Plot_Generation_ver2.Rmd", package="icreport")

    outFile = paste(output.path,"/",prefix,"_PCA_summary.html",sep="")

    rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document")

}
