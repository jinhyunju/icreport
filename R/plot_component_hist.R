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
#' @export
#'
#' @examples
#' R code here showing how your function works
plot_component_hist <- function(s,plot.title){
  G <- length(s)
  plot.sig <- data.frame(idx = c(1:G), sig = s)
  p <- ggplot(plot.sig, aes(x = sig)) +geom_histogram(aes(y = ..density..)) + xlab("Gene Weights") +
    labs(title = plot.title) + geom_density() + ylab("")
  return(p)

}
