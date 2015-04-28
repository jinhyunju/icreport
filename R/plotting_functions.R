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
#' @examples
#' R code here showing how your function works
#' @export
plot_component_hist <- function(s, plot.title){
  G <- length(s)
  plot.sig <- data.frame(idx = c(1:G), sig = s)
  p <- ggplot(plot.sig, aes(x = sig)) +geom_histogram(aes(y = ..density..)) + xlab("Gene Weights") +
    labs(title = plot.title) + geom_density() + ylab("") + coord_flip()
  return(p)

}


#' Plotting the gene weights of a single component
#'
#' ggplot function for plotting a single component
#'
#' @param s A single component  \code{s}
#' @param plot.title Character title of the plot\code{plot.title}
#' @return ggplot object
#'
#' @keywords keywords
#'
#'
#' @examples
#' R code here showing how your function works
#' @export
plot_component <- function(s,
                           plot.title,
                           peaks = FALSE,
                           peakresult = NULL,
                           gene.names = NULL){
  if(peaks == FALSE){
    G <- length(s)
    #plot.sig <- data.frame(idx = seq(1,G,1), sig = s)
    plot.sig <- data.frame(idx = rank(-s, ties.method = "first"), sig = s)
    p <- ggplot(plot.sig, aes(x=idx, y= sig)) +
      geom_linerange(aes(ymin=0, ymax=sig)) +
      xlab("Gene index") +
      ylab("Gene Weights")+ labs(title = plot.title)
    return(p)
  } else {
    G <- length(s)
    colnames(s) <- "sig"
    significant.peaks <- names(peakresult)
    n.peaks <- length(peakresult)
    peak.idx <- 1 * (gene.names %in% significant.peaks)
    #plot.sig <- data.frame(idx = c(1:G), sig = s, peaks = peak.idx)
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
#'
#' @examples
#' R code here showing how your function works
plot_coeff_w_legend <- function(info_input, k.col){

  if(k.col != "none"){
    ggplot(info_input, aes(x=idx, y= IC)) +
      geom_point(aes_string(colour = k.col),shape=20, size = 3) +
      xlab("sample") + ylab("IC coefficient")
  } else {
    ggplot(info_input, aes(x=idx, y= IC)) +
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
#' @export
#'
#' @examples
#' R code here showing how your function works
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

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
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
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
#' @export
plot_component_chr <- function(s,
                               geneinfo.df,
                               x.axis,
                               plot.title,
                               peaks = FALSE,
                               peakresult = NULL){
  geneinfo.df$loading <- s[as.character(geneinfo.df$phenotype),]
  geneinfo.df$peaks <- 1 * ( as.character(geneinfo.df$phenotype) %in% names(peakresult) )
  #1 * (abs(geneinfo.df$loading) > (2 * sd(geneinfo.df$loading)))
  n.peaks <- sum(geneinfo.df$peaks)
  plot.title <- paste(plot.title,"_",n.peaks,"peaks", sep = " ")

#  rect_left <- x.axis[1,seq(1, ncol(x.axis), by = 2)]
#  rect_right <- x.axis[2,seq(1, ncol(x.axis), by = 2)]
  rect_left <- x.axis[1,]
  rect_right <- x.axis[2,]

  rects1 <- data.frame(
    xstart = rect_left,
    xend = rect_right,
    idx = rep(0),
    loading = rep(0),
    fill.col = factor(rep(c(0,1), length.out = length(rect_left)))
  )


  p <- ggplot(geneinfo.df, aes(x = idx, y = loading)) +
    geom_rect(data = rects1, aes(xmin = xstart,
                                 xmax = xend,
                                 ymin = -Inf,
                                 ymax = Inf,
                                 fill = fill.col), alpha = 0.35) +
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
