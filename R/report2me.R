#' Generating ICA or PCA summary reports.
#'
#' Generating a HTML report from a ICA or PCA list object.
#'
#' @param input ICA or PCA result list created by either \code{gene_expr_ica()} or \code{gene_expr_pca()}.
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
#' @import rmarkdown
#' @import gtools
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
report2me <- function(input = NULL,
                      n.comps = NULL,
                      prefix = NULL,
                      geneinfo.df = NULL,
                      output.path = NULL, file.ext = "png"){

    if(!exists("input")){
        stop("Please specify the input to generate a report. \n")
    }

    if(!exists("prefix")){
        stop("Please specify the prefix of the output file \n")
    }

    if(is.null(output.path)){
        message("Output path is not specified, using current working directory \n")
        output.path <- getwd()

    }

    method <- attr(input, 'method')

    if(!is.null(geneinfo.df)){
      geneinfo.df$pheno_chr <- as.character(geneinfo.df$pheno_chr)
      chromosomes <- unique(geneinfo.df$pheno_chr)
      chr.ordered <- gtools::mixedsort(chromosomes)
      geneinfo.df <- geneinfo.df[,c("phenotype","pheno_chr","pheno_start","pheno_end")]
      geneinfo.df <- geneinfo.df[!duplicated(geneinfo.df),]
    }

    if( method == "ica"){

        ica.result <- input
        info.df <- ica.result$info.df
        data.set <- prefix

        rm(input, prefix)

        if(!is.null(geneinfo.df)){
          message("Checking if gene information exists for every gene\n")

          geneposition.df <- subset(geneinfo.df, phenotype %in% rownames(ica.result$S))

          n.genesinfo <- sum(unique(geneposition.df$phenotype) %in% rownames(ica.result$S))

          if(n.genesinfo != nrow(ica.result$S)){

            message("Missing some gene position information, NA category created")
            gene.names <- data.frame("phenotype" = rownames(ica.result$S))
            merged.geneinfo <- merge(gene.names, geneposition.df, all = TRUE)

            # order entries according to chr number and position
            merged.geneinfo <- merged.geneinfo[order(match(merged.geneinfo$pheno_chr, chr.ordered), merged.geneinfo$pheno_start),]

            #removing potential duplicate entries
            geneposition.df <- merged.geneinfo[!duplicated(merged.geneinfo$phenotype),]
          } else {

            message("All genes have chromosome and position information\n")
            geneposition.df <- geneposition.df[order(match(geneposition.df$pheno_chr, chr.ordered), geneposition.df$pheno_start),]
            geneposition.df <- geneposition.df[!duplicated(geneposition.df$phenotype),]
          }


          ica.result$ordered.geneinfo <- geneposition.df
          ica.result$ordered.geneinfo$idx <- c(1:nrow(ica.result$ordered.geneinfo))
          rm(geneposition.df)

        }




        markdown.file <- system.file("templates/Component_Visualization_Report.Rmd", package="icreport")

        outFile = paste(output.path,"/",data.set,"_ICA_summary.html",sep="")

        suppressMessages(rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document"))

    } else if (method == "pca"){

        pca.result <- input

        if(is.null(n.comps)){
          plot.pc <- dim(pca.result$x)[2]
        } else {
          plot.pc <- n.comps
        }

        if(!is.null(geneinfo.df)){
          cat("Checking if gene information exists for every gene\n")

          geneposition.df <- subset(geneinfo.df, phenotype %in% rownames(pca.result$rotation))

          pca.result$ordered.geneinfo <- geneposition.df[order(nchar(geneposition.df$pheno_chr),geneposition.df$pheno_chr),]
          pca.result$ordered.geneinfo$idx <- c(1:nrow(pca.result$ordered.geneinfo))

          if(nrow(pca.result$ordered.geneinfo) != nrow(pca.result$rotation)){
            stop("Missing gene position information");
          } else {
            cat("All genes have chromosome and position information\n")
          }

        }


        data.set <- prefix
        #markdown.file <- system.file("templates/PCA_Plot_Generation_ver2.Rmd", package="icreport")
        markdown.file <- system.file("templates/Component_Visualization_Report.Rmd", package="icreport")
        outFile = paste(output.path,"/",prefix,"_PCA_summary.html",sep="")

        rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document")


    }


}
