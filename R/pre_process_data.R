#' Pre-processing of the data
#'
#' The function takes in a g x N matrix, where g is the number of genes and N is the number of samples,
#' and removes any 0 variance entries followed by mean centering.
#' If scaling option is true it will scale all rows to have unit standard deviation.
#'
#' @param phenotype.mx Phenotype matrix with a dimension of g X N  \code{phenotype.mx}
#' @return Pre-processed phenotype matrix
#' @keywords keywords
#'
#'
#' @examples
#' raw.data <- rbind(rnorm(100,2,3),rnorm(100,10,2))
#' centered <- pre_process_data(raw.data, scale.pheno = FALSE)
#' apply(centered, 1, mean)
#' apply(centered, 1, sd)
#' centered.scaled <- pre_process_data(raw.data, scale.pheno = TRUE)
#' apply(centered.scaled, 1, mean)
#' apply(centered.scaled, 1, sd)
#' @export
pre_process_data <- function(phenotype.mx, scale.pheno = TRUE){
    var.vec <- apply(phenotype.mx, 1, var)
    var0.count <- sum(var.vec == 0)


    # Remove 0 variance genes if existing
    if(var0.count != 0){
        message("- Removing ",var0.count," 0 variance genes \n")
        phenotype.mx <- phenotype.mx[-c(which(var.vec == 0)),]
    } else {
        message("- No 0 variance genes in dataset \n")
    }


    if(scale.pheno == TRUE){
        message("- Centering and Scaling Phenotype Matrix \n")
        # scale phenotype matrix to have 0 mean and 1 sd
        phenotype.mx <- t(apply(phenotype.mx, 1, function(x) (x - mean(x))/sd(x)))
    } else {
        message("- Centering Phenotype Matrix \n")
        phenotype.mx <- t(apply(phenotype.mx, 1, function(x) (x - mean(x))))
    }

    message("- Dimension of Phenotype Matrix = [",dim(phenotype.mx)[1],",",dim(phenotype.mx)[2],"]\n")
    return(phenotype.mx)
}

