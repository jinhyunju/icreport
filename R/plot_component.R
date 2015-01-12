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
#' @export
#'
#' @examples
#' R code here showing how your function works
plot_component <- function(s,plot.title,peaks = FALSE, peakresult = NULL, gene.names = NULL){
  if(peaks == FALSE){
    G <- length(s)
    plot.sig <- data.frame(idx = seq(1,G,1), sig = s)
    p <- ggplot(plot.sig, aes(x=idx, y= sig)) +
      geom_linerange(aes(ymin=0, ymax=sig)) +
      xlab("Gene index") +
      ylab("Gene Weights")+ labs(title = plot.title)
    return(p)
  } else {
    G <- length(s)
    significant.peaks <- names(peakresult)
    n.peaks <- length(peakresult)
    peak.idx <- 1 * (gene.names %in% significant.peaks)
    plot.sig <- data.frame(idx = c(1:G), sig = s, peaks = peak.idx)
    plot.title <- paste(plot.title,"_",n.peaks,"peaks", sep = " ")
    p <- ggplot(plot.sig, aes(x=idx, y= sig)) + geom_linerange(aes(ymin=0, ymax=sig,colour = factor(peaks))) +
      scale_y_continuous(expand=c(0,0))+scale_color_manual(values=c("black", "red")) + xlab("Gene index") + ylab("Gene Weights")+
      labs(title = plot.title) +theme(legend.position="none")
    return(p)
  }

}
