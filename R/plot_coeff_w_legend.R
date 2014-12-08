#' Plotting the coefficients for a single IC
#'
#' ggplot function for plotting coefficients for a single IC
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
plot_coeff_w_legend <- function(info_input,k.col,plot.title){

    if(k.col != "none"){
        ggplot(info_input, aes(x=idx, y= IC)) +
        geom_point(aes_string(colour = k.col),shape=20, size = 3) +
        xlab("sample") + ylab("IC coefficient") +
        labs(title = plot.title)
    } else {
        ggplot(info_input, aes(x=idx, y= IC)) +
        geom_point(aes(color = "black"),shape=20, size = 3) +
        xlab("sample") + ylab("IC coefficient") +
        labs(title = plot.title)
    }

}
