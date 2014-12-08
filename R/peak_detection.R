#' Labeling peaks for each independent component
#'
#' The function takes in the A matrix from the ica.result object
#' with a dataframe of measured covariates to test the association
#' between them.
#'
#' @param s A single column of the S matrix with dimensions 1 x g  \code{s}
#' @return list with entries "peaks" and "N". "peaks" has the position of peak genes
#' in the decreasing order of magnitude. "N" has the total number of peak genes.
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
peak_detection <- function(ica.input.s){
  peaks.idx <- which(abs(ica.input.s) > (2 * sd(ica.input.s)))
  peaks <- sort(abs(ica.input.s[peaks.idx]), decreasing = T)
  return(list("peaks" = peaks, "N"= length(peaks)))
}
