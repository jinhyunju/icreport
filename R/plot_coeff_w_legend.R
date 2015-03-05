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
plot_coeff_w_legend <- function(info_input,k.col){

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
