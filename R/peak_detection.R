#' Labeling peaks for each independent component
#'
#' The function takes in the A matrix from the ica.result object
#' with a dataframe of measured covariates to test the association
#' between them.
#'
#' @param s A single column of the S matrix with dimensions 1 x g  \code{s}
#' @return A vector that contains the gene weights of the signficant peaks for a single independent component sorted
#' in the decreasing order of absolute magnitude.
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
peak_detection <- function(ica.input.s){
  peaks.idx <- which(abs(ica.input.s) > (2 * sd(ica.input.s))) # get the peak indexes
  peaks <- ica.input.s[names(sort(abs(ica.input.s[peaks.idx]), decreasing = T))] # sort them in decreasing order of absolute magnitude
  return(peaks)
}
