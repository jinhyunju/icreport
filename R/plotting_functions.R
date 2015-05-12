#' Plotting a histogram of gene weights of a single component
#'
#' ggplot function for plotting a histogram for a single component
#'
#' @param s A single component  \code{s}
#' @param plot.title Character title of the plot\code{plot.title}
#' @return ggplot object
#'
#' @keywords keywords
#'
#'
#' @export
plot_component_hist <- function(s, plot.title){
  G <- length(s)
  plot.sig <- data.frame(idx = c(1:G), sig = s)
  p <- ggplot(plot.sig, aes_string(x = "sig")) +geom_histogram(aes_string(y = "..density..")) + xlab("Gene Weights") +
    labs(title = plot.title) + geom_density() + ylab("") + coord_flip()
  return(p)

}


#' Plotting the gene weights of a single component
#'
#' ggplot function for plotting a single component
#'
#' @param s A single component
#' @param plot.title Character title of the plot
#' @param peaks TRUE shows peaks in red, FALSE leaves everything black
#' @param peakresult Pre-calculated peak positions, if not supplied peaks will be
#'        defined as values lying out of 2 standard deviations.
#' @return ggplot object
#'
#' @keywords keywords
#'
#' @export
plot_component <- function(s,
                           plot.title,
                           peaks = FALSE,
                           peakresult = NULL){
  sig <- NULL
  idx <- NULL # both to avoid R CMD check NOTES
  if(peaks == FALSE){
    G <- length(s)

    plot.sig <- data.frame(idx = rank(-s, ties.method = "first"), sig = s)
    p <- ggplot(plot.sig, aes_string(x= "idx", y= "sig")) +
      geom_linerange(aes(ymin=0, ymax=sig)) +
      xlab("Gene index") +
      ylab("Gene Weights")+ labs(title = plot.title)
    return(p)
  } else {
    G <- length(s)

    colnames(s) <- "sig"

    significant.peaks <- names(peakresult)

    n.peaks <- length(peakresult)

    peak.idx <- 1 * (rownames(s) %in% significant.peaks)

    plot.sig <- data.frame(idx = rank(-s, ties.method = "first"), sig = s, peaks = peak.idx)
    plot.title <- paste(plot.title,"_",n.peaks,"peaks", sep = " ")
    p <- ggplot(plot.sig, aes(x=idx, y= sig)) +
      geom_linerange(aes(ymin=0, ymax=sig, colour = factor(peaks))) +
      scale_y_continuous(expand=c(0,0))+scale_color_manual(values=c("black", "red")) +
      xlab("Genes sorted by weights") + ylab("Gene Weights")+
      labs(title = plot.title) +theme(legend.position="none")
    return(p)
  }

}

#' Plotting the coefficients for a single IC
#'
#' ggplot function for plotting coefficients for a single IC
#'
#' @param info_input A dataframe that has the index and values of ICs with covariate information saved.
#' @param k.col The column name to use for plot coloring.
#' @return ggplot object
#'
#' @keywords keywords
#'
#' @export
plot_coeff_w_legend <- function(info_input, k.col){

  if(k.col != "none"){
    plot.factor <- factor(info_input[,k.col])
    number.of.levels <- length(levels(plot.factor))
      if(number.of.levels < 10){
        ggplot(info_input, aes_string(x="idx", y= "IC")) +
          geom_point(aes_string(colour = k.col),shape=20, size = 3) +
          xlab("sample") + ylab("IC coefficient")
      } else {
        ggplot(info_input, aes_string(x= "idx", y= "IC")) +
          geom_point(aes_string(colour = k.col),shape=20, size = 3) +
          xlab("sample") + ylab("IC coefficient") + theme(legend.position = "none")
      }
    } else {
    ggplot(info_input, aes_string(x="idx", y= "IC")) +
      geom_point(aes(color = "black"),shape=20, size = 3) +
      xlab("sample") + ylab("IC coefficient")
  }

}

#' Plotting multiple ggplots
#'
#' Equivalent function for ggplots to setting par(mfrow=c(x,y)) for
#' basic R plots.
#' Adapted from : http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#'
#' @param plotlist list of plots
#' @param cols number of columns
#' @param layout layout matrix for custom layout
#' @return single plot object
#'
#' @keywords keywords
#'
#' @import grid
#'
#' @export
multiplot <- function(plotlist=NULL, layout=NULL, cols=1) {

  # Make a list from the ... arguments and plotlist
  #plots <- c(list(...), plotlist)
  plots <- plotlist
  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Adding customized black and white theme to ggplot objects and controls the font size.
#'
#' The function takes a ggplot object as an input and modifies the theme of it
#'
#' @param inputplot A ggplot object.  \code{s}
#' @return A ggplot object with a custom theme added
#' @export
ggplot_add_theme <- function(inputplot){
  inputplot <- inputplot + theme(axis.text = element_text(size=12),
                                 axis.title=element_text(size=14,face="bold"),
                                 title = element_text(size=16,face="bold"),
                                 legend.title = element_text(size= 10, face="bold"),
                                 legend.text = element_text(size = 8),
                                 legend.key = element_rect(fill = NA,color = "white"),
                                 panel.grid.major = element_line(colour = "grey80", linetype="dashed"),
                                 panel.background = element_rect(fill = NA,color = "grey80"))
  return(inputplot)
}

#' Plotting function that indicates gene position on chromosome
#'
#' The function generates a plot with gene weights for individual ICs according
#' to their position on each chromosome.
#'
#' @param s A single component
#' @param geneinfo.df Gene position information dataframe
#' @param x.axis Positions of the chromosome borders for coloring.
#' @param plot.title Character title of the plot
#' @param peaks TRUE shows peaks in red, FALSE leaves everything black
#' @param peakresult Pre-calculated peak positions, if not supplied peaks will be
#'        defined as values lying out of 2 standard deviations.
#'
#' @export
plot_component_chr <- function(s,
                               geneinfo.df,
                               x.axis,
                               plot.title,
                               peaks = TRUE,
                               peakresult = NULL){
  geneinfo.df$loading <- s[as.character(geneinfo.df$phenotype),]
  loading <- NULL # avoid R CMD check NOTES

  if(peaks == TRUE){
    if(is.null(peakresult)){
      geneinfo.df$peaks <- 1 * (abs(geneinfo.df$loading) > (2 * sd(geneinfo.df$loading)))
    } else {
      geneinfo.df$peaks <- 1 * ( as.character(geneinfo.df$phenotype) %in% names(peakresult) )
    }

  } else {
      geneinfo.df$peaks <- rep(0, nrow(geneinfo.df))
  }


  n.peaks <- sum(geneinfo.df$peaks)
  plot.title <- paste(plot.title,"_",n.peaks,"peaks", sep = " ")

  rect_left <- x.axis[1,]
  rect_right <- x.axis[2,]

  rects1 <- data.frame(
    xstart = rect_left,
    xend = rect_right,
    idx = rep(0),
    loading = rep(0),
    fill.col = factor(rep(c(0,1), length.out = length(rect_left)))
  )


  p <- ggplot(geneinfo.df, aes_string(x = "idx", y = "loading")) +
    geom_rect(data = rects1, aes_string(xmin = "xstart",
                                 xmax = "xend",
                                 ymin = "-Inf",
                                 ymax = "Inf",
                                 fill = "fill.col"), alpha = 0.35) +
    geom_linerange(aes(ymin = 0, ymax = loading, colour = factor(peaks))) +
    theme_bw() +
    #geom_vline(xintercept=c(0,x.axis[2,]), linetype="dotted") +
    scale_x_continuous("Gene Chromosome Location",breaks=x.axis[3,],labels=unique(geneinfo.df$pheno_chr))+
    scale_y_continuous(expand = c(0, 0)) + labs(title = plot.title) +
    scale_color_manual(values=c("grey", "red")) +
    scale_fill_manual(values=c("skyblue", NA)) +
    ylab("Gene Weights") +
    theme(legend.position="none")

  return(p)
}

#' Function that generates the x axis ticks for chromosome positions
#'
#' The function takes in a gene position info data frame with the following columns:
#' "phenotypes" holding the names of the genes used in the analysis (rownames of phenotype.mx)
#' "pheno_chr" holding the chromosome each gene belongs to
#' "pheno_start" the starting position of a given gene
#' "pheno_end" the end position of a given gene
#' "idx" a non-overlapping index that shows orders the genes based on their position, which goes from
#' 1 to number of genes.
#'
#' @param geneinfo.df dataframe that holds the position of the genes.
#'
#' @export
chr_axis_creator <- function(geneinfo.df){
  if(is.null(geneinfo.df$idx)){
    stop("index is missing in the geneinfo dataframe. Please assign indexes first")
  }
  pheno_chr <- NULL # to get rid of R CMD check NOTES
  chromosomes <- unique(geneinfo.df$pheno_chr)
  x.axis <- matrix(0, nrow = 3, ncol = length(chromosomes))

  for(k in 1:length(chromosomes)){
    j <- chromosomes[k]

    if(is.na(j)){
      x.axis[1,k] <- max(x.axis[1,]) + 1
      x.axis[2,k] <- max(geneinfo.df[,"idx"])
    } else {
      x.axis[1,k] <- min(subset(geneinfo.df, pheno_chr == j)[,"idx"])
      x.axis[2,k] <- max(subset(geneinfo.df, pheno_chr == j)[,"idx"])
    }


  }

  x.axis[3,] <- (x.axis[1,] + x.axis[2,]) / 2
  return (x.axis)
}
