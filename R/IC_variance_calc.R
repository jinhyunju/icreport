#' Calculating the variation explained by each component
#'
#' The function takes in a single row of the A matrix and a
#' single column of the S matrix and calculates the sums of squares.
#'
#' @param s A single column of the S matrix with dimensions 1 x g  \code{s}
#' @param a A single row of the A matrix with dimensions N x 1 \code{a}
#' @return A single number of how much variantion is explained by a single IC.
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
IC_variance_calc <- function(s, a){
  var.IC <- sum(  (as.matrix(s) %*% t(as.matrix(a) ) )^2)
  return(var.IC)
}
