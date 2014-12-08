#' Automated ICA report generator
#'
#' Generating a HTML report from a ICA list object.
#'
#' @param ica.result ICA result list created by the function run_ica.   \code{ica.result}
#' @param prefix Output filename prefix. The output file will be named
#'        "prefix_ICA_summary.html". \code{prefix}
#' @param output.path Directory path for generating the output HTML file.
#'        default is set to current working directory. \code{output.path}
#' @return output HTML report that
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
ic_report <- function(ica.result,prefix,output.path = NULL, file.ext = "png"){
    if(!exists("prefix")){
        cat("Please specify the prefix of the output file \n")
        break;
    }

    if(is.null(output.path)){
        cat("Output path is not specified, using current working directory \n")
        output.path <- getwd()

    }

    info.df <- ica.result$info.df

    data.set <- prefix
    markdown.file <- system.file("templates/ICA_Plot_Generation_ver2.Rmd", package="icreport")

    outFile = paste(output.path,"/",prefix,"_ICA_summary.html",sep="")

    rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document")

}
