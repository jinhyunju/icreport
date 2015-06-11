#' Testing the association between independent component coefficients and known covariates.
#'
#' The function takes in the A matrix from the ica.result object
#' with a dataframe of measured covariates to test the association
#' between them.
#'
#' @param input.A The mixing matrix estimated through ICA
#' @param info.input A dataframe that holds the measured covariates for each sample
#' @param covar.names names of the covariates that should be used for association testing
#' @return output A matrix holding the p-values for each indepenedent component and covariate pair.
#' @keywords keywords
#'
#' @export
#'
component_association_test <- function(input.A, info.input, covar.names){

  n.components <- dim(input.A)[1]
  n.covariates <- length(covar.names)

  pval.mx <- matrix(NA, nrow = n.components, ncol = n.covariates)
  colnames(pval.mx) <- covar.names

  covar.df <- data.frame(info.input[colnames(input.A),covar.names])

  for( i in 1: dim(input.A)[1]){
    for (j in 1:length(covar.names)){
      analysis.df <- data.frame("IC"=input.A[i,],"Covar" = covar.df[,j])
      anova.fit <- aov(IC ~ Covar, data=analysis.df)
      model.fstat <- summary.lm(anova.fit)$fstatistic
      p.val <- pf(model.fstat[1],model.fstat[2],model.fstat[3], lower.tail = F, log.p = F)
      pval.mx[i,j] <- p.val
    }
  }
  return(pval.mx)
}

#' Labeling peaks for each independent component
#'
#' The function takes in the A matrix from the ica.result object
#' with a dataframe of measured covariates to test the association
#' between them.
#'
#' @param ica.input.s A single column of the S matrix with dimensions 1 x g  \code{s}
#' @return A vector that contains the gene weights of the signficant peaks for a single independent component sorted
#' in the decreasing order of absolute magnitude.
#' @keywords keywords
#'
#' @export
peak_detection <- function(ica.input.s){
  peaks.idx <- which(abs(ica.input.s) > (2 * sd(ica.input.s))) # get the peak indexes
  peaks <- ica.input.s[names(sort(abs(ica.input.s[peaks.idx]), decreasing = T))]
  # sort them in decreasing order of absolute magnitude
  return(peaks)
}

#' Calculating the variation explained by each component
#'
#' The function takes in a single row of the ICA result A matrix and a
#' single column of the ICA result S matrix and calculates the sums of squares.
#'
#' @param s A single column of the S matrix with dimensions 1 x g  \code{s}
#' @param a A single row of the A matrix with dimensions N x 1 \code{a}
#' @return A single number of how much variantion is explained by a single IC.
#'
#' @export
IC_variance_calc <- function(s, a){
  var.IC <- sum(  (as.matrix(s) %*% t(as.matrix(a) ) )^2)
  return(var.IC)
}
