#' Testing the association between independent component coefficients and known covariates.
#'
#' The function takes in the A matrix from the ica.result object
#' with a dataframe of measured covariates to test the association
#' between them.
#'
#' @param A The mixing matrix estimated through ICA   \code{A}
#' @param info.df A dataframe that holds the measured covariates for each sample \code{info.df}
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

      covar.mx <- as.matrix(info.input[colnames(input.A),covar.names])

      for( i in 1: dim(input.A)[1]){
        for (j in 1:length(covar.names)){
          analysis.df <- data.frame("IC"=input.A[i,],"Covar" = covar.mx[,j])
          anova.fit <- aov(IC ~ Covar, data=analysis.df)
          model.fstat <- summary.lm(anova.fit)$fstatistic
          p.val <- pf(model.fstat[1],model.fstat[2],model.fstat[3], lower.tail = F, log.p = F)
          pval.mx[i,j] <- p.val
        }
      }
      return(pval.mx)
}
