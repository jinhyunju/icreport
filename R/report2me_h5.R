#' Generating ICA or PCA summary reports.
#'
#' Generating a HTML report from a ICA or PCA list object.
#'
#' @param input ICA or PCA result list created by either \code{gene_expr_ica()} or \code{gene_expr_pca()}.
#' @param n.comps Number of principal components to plot.
#'        Default is set to plot every PC. Only used for PCA plotting not for ICA plotting.
#' @param prefix Output filename prefix. The output file will be named
#'        "prefix_ICA_summary.html".
#' @param geneinfo.df Dataframe that contains positions of the genes. Column names should be
#'        "pheno_chr" for chromosomes, "pheno_start" for starting position and "pheno_end" for
#'        ending positions.
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
#' @import rhdf5
#' @import markdown
#'
#' @export
h5report <- function(input_h5 = NULL,
                     output.path = NULL){

    if(!exists("input_h5")){
        stop("Missing hdf5 file input for report. \n")
    }

    if(is.null(output.path)){
        message("Output path is not specified, using current working directory \n")
        output.path <- getwd()

    }

    gene_ids <- h5read(input_h5, "phenotypes/col_info/id")
    gene_chr <- h5read(input_h5, "phenotypes/col_info/pheno_chr")

    snp_ids <- h5read(input_h5, "genotypes/col_info/id")
    snp_chr <- h5read(input_h5, "genotypes/col_info/geno_chr")
    snp_position <- h5read(input_h5, "genotypes/col_info/geno_pos")


    covars <- h5read(input_h5, "covars/matrix")
    colnames(covars) <- h5read(input_h5, "covars/col_info/id")
    sample_ids <- h5read(input_h5, "phenotypes/row_info/id")

    message("Running PCA on phenotypes")

    phenotypes <- h5read(input_h5, "phenotypes/matrix")
    colnames(phenotypes) <- gene_ids
    rownames(phenotypes) <- sample_ids

    pheno_pca <- gene_expr_pca(t(phenotypes))

    message("Running PCA on genotypes")
    genotypes <- h5read(input_h5, "genotypes/matrix")
    rownames(genotypes) <- sample_ids
    colnames(genotypes) <- snp_ids


    geno_pca <- gene_expr_pca(t(genotypes))

    n_samples <- nrow(phenotypes)
    n_pheno <- ncol(phenotypes)
    n_geno <- ncol(genotypes)


    AF <- colMeans(genotypes) / 2
    MAF <- pmin(AF, 1 -AF)

    snp_info <- data.frame("id" = snp_ids,
                           "chr" = factor(snp_chr, levels = gtools::mixedsort(snp_chr)),
                           "MAF" = MAF,
                           "position" = snp_position)
    message("Calculating gene variance")
    gene_var <- apply(phenotypes, 2, var)

    gene_info <- data.frame("id" = gene_ids,
                            "chr" = factor(gene_chr, levels = gtools::mixedsort(gene_chr)),
                            "var" = gene_var)


    data.set <- gsub("[.]h5", "", input_h5)
    markdown.file <- system.file("templates/h5_report.Rmd", package="icreport")
    #markdown.file <- system.file("./Component_Visualization_Report.Rmd", package="icreport")
    outFile = paste(output.path,"/",data.set,"_Data_Report.html",sep="")

    out <- knitr::knit(markdown.file)
    markdown::markdownToHTML(out, outFile)
}

#' Simplified Plotting function that indicates gene position on chromosome
#'
#' The function generates a plot with gene weights for individual ICs according
#' to their position on each chromosome.
#'
#' @param s A single component
#' @param geneinfo_df Gene position information dataframe with columns "chr"
#' @param plot.title Character title of the plot
#'
#' @export

plot_component_grid <- function(s,geneinfo_df,plot.title){

  
  geneinfo_df$loading <- s
  geneinfo_df$peaks <- 1 * (abs(geneinfo_df$loading) > (2 * sd(geneinfo_df$loading)))
#  geneinfo_df <- with(geneinfo_df, geneinfo_df[order(chr, loading),])
  geneinfo_df$idx <- c(1:nrow(geneinfo_df))
  chromosomes <- unique(geneinfo_df$chr)
  
  x.axis <- matrix(0, nrow = 3, ncol = length(chromosomes))
  
  for(k in 1:length(chromosomes)){
    j <- chromosomes[k]
    
    if(is.na(j)){
      x.axis[1,k] <- max(x.axis[1,]) + 1
      x.axis[2,k] <- max(geneinfo_df[,"idx"])
    } else {
      x.axis[1,k] <- min(subset(geneinfo_df, chr == j)[,"idx"])
      x.axis[2,k] <- max(subset(geneinfo_df, chr == j)[,"idx"])
    }
    
    
  }
  
  x.axis[3,] <- (x.axis[1,] + x.axis[2,]) / 2
  
  plot.title <- paste0(plot.title,"_",sum(geneinfo_df$peaks),"peaks")
  
  rect_left <- x.axis[1,]
  rect_right <- x.axis[2,]
  
  rects1 <- data.frame(
    xstart = rect_left,
    xend = rect_right,
    idx = rep(0),
    loading = rep(0),
    fill.col = factor(rep(c(0,1), length.out = length(rect_left)))
  )
  
  
  p <- ggplot(geneinfo_df, aes_string(x = "idx", y = "loading")) + 
              geom_linerange(aes(ymin = 0, ymax = loading, colour = factor(peaks))) + 
 #             geom_rect(data = rects1, aes_string(xmin = "xstart",
#                                                  xmax = "xend",
#                                                  ymin = "-Inf",
#                                                  ymax = "Inf",
#                                                  fill = "fill.col"), alpha = 0.35) + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ 
    facet_grid(.~chr, scales = "free_x", space = "free_x") + labs(title = plot.title) +
#    scale_x_continuous("Gene Chromosome Location",breaks=x.axis[3,],labels=unique(geneinfo_df$chr))+
#    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(values=c("grey", "red")) +
    scale_fill_manual(values=c("skyblue", NA)) +
    ylab("Gene Weights") + xlab("") +
    theme(legend.position="none")
  
  return(p)
}

