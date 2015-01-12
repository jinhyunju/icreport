#' Generating ICA or PCA summary reports.
#'
#' Generating a HTML report from a ICA or PCA list object.
#'
#' @param input ICA or PCA result list created by either \code{gene_expr_ica()} or \code{gene_expr_pca()}.
#' @param mode Specifying whether the object is a "ica" object or "pca" object. Default is set to "ica".
#' @param n.comps Number of principal components to plot. Default is set to plot every PC. Only used for PCA plotting not for ICA plotting.
#' @param prefix Output filename prefix. The output file will be named
#'        "prefix_ICA_summary.html".
#' @param output.path Directory path for generating the output HTML file.
#'        default is set to current working directory.
#' @param file.ext File extension to be used for saved plots. Default is set to png for html reports.
#'        Note that if you use pdf plots for html files they will not show up in your report.
#' @return output HTML report.
#' @keywords keywords
#'
#' @import ggplot2
#' @import knitr
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
report2me <- function(input = NULL, mode = "ica", n.comps = NULL, prefix = NULL, output.path = NULL, file.ext = "png"){

    if(!exists("input")){
      cat("Please specify the input to generate a report. \n")
      break;
    }

    if(!exists("prefix")){
        cat("Please specify the prefix of the output file \n")
        break;
    }

    if(is.null(output.path)){
        cat("Output path is not specified, using current working directory \n")
        output.path <- getwd()

    }


    if( mode == "ica"){

        ica.result <- input
        info.df <- ica.result$info.df
        data.set <- prefix

        markdown.file <- system.file("templates/ICA_Plot_Generation_extended.Rmd", package="icreport")

        outFile = paste(output.path,"/",prefix,"_ICA_summary.html",sep="")

        suppressMessages(rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document"))

    } else if (how == "pca"){

        pca.result <- input

        if(is.null(n.comps)){
          plot.pc <- dim(pca.result$x)[2]
        } else {
          plot.pc <- n.comps
        }

        data.set <- prefix
        markdown.file <- system.file("templates/PCA_Plot_Generation_ver2.Rmd", package="icreport")

        outFile = paste(output.path,"/",prefix,"_PCA_summary.html",sep="")

        rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document")


    }


}