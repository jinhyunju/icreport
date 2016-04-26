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

    markdown.file <- system.file("templates/ICA_Visualization.Rmd", package="icreport")

    outFile = paste(output.path,"/",prefix,"_ICA_summary.html",sep="")

    suppressMessages(rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document"))



}
