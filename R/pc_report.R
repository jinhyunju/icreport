#' Automated PCA report generator
#'
#' Generating a HTML report from a PCA object created with the function prcomp().
#'
#' @param pca.result PCA result list created by the function prcomp().   \code{pca.result}
#' @param plot.pc Number of principal components to plot. Default is set to plot every PC. \code{plot.pc}
#' @param prefix Output filename prefix. The output file will be named
#'        "prefix_ICA_summary.html". \code{prefix}
#' @param output.path Directory path for generating the output HTML file.
#'        default is set to current working directory. \code{output.path}
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
pc_report <- function(pca.result = NULL, plot.pc = NULL, prefix = NULL, output.path = NULL, file.ext = "png"){

    if(is.null("pca.result")){
      cat("pca.result object is missing. Please generate a pca.result object using the function gene_expr_pca() \n")
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

    if(is.null(plot.pc)){
        plot.pc <- dim(pca.result$x)[2]
    }

    data.set <- prefix
    markdown.file <- system.file("templates/PCA_Plot_Generation_ver2.Rmd", package="icreport")

    outFile = paste(output.path,"/",prefix,"_PCA_summary.html",sep="")

    rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document")

}
