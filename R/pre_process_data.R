#' Pre-processing of the data
#'
#' The function takes in the phenotype matrix
#' and removes any 0 variance entries. Also it mean centers the data.
#' If scaling option is true it will scale all the rows to have unit standard deviation.
#'
#' @param phenotype.mx Phenotype matrix with a dimension of g X N  \code{phenotype.mx}
#' @return Pre-processed phenotype matrix
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
pre_process_data <- function(phenotype.mx, scale.pheno = TRUE){
    var.vec <- apply(phenotype.mx, 1, var)
    var0.count <- sum(var.vec == 0)


    # Remove 0 variance genes if existing
    if(var0.count != 0){
        cat("Removing ",var0.count," 0 variance genes \n")
        phenotype.mx <- phenotype.mx[-c(which(var.vec == 0)),]
    } else {
        cat("No 0 variance genes in dataset \n")
    }


    if(scale.pheno == TRUE){
        cat("Centering and Scaling Phenotype Matrix \n")
        # scale phenotype matrix to have 0 mean and 1 sd
        phenotype.mx <- t(apply(phenotype.mx, 1, function(x) (x - mean(x))/sd(x)))
    } else {
        cat("Centering Phenotype Matrix \n")
        phenotype.mx <- t(apply(phenotype.mx, 1, function(x) (x - mean(x))))
    }

    cat("Dimension of Phenotype Matrix = [",dim(phenotype.mx)[1],dim(phenotype.mx)[2],"]\n")
    return(phenotype.mx)
}

