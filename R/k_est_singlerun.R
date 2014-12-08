#' Estimating the number of components (k)
#'
#' Estimating the number of components to be estimated
#' through the replication ratio between two independent runs of ICA
#'
#' @param phenotype.mx The matrix that holds the gene expression profile   \code{phenotype.mx}
#' @param max.iter Number of maximum iterations, default is set to 10 \code{max.iter}
#' @param run index \code{run.index}
#' @return output A matrix holding the p-values for each indepenedent component and covariate pair.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
k_est_singlerun <- function(phenotype.mx = NULL,max.iter = 10,run.index = 1){
  summary.vec <- NULL
  k.vec <- NULL
  if( is.null(phenotype.mx)){
    cat("Phenotype matrix is missing \n")
    stop;
  }
  # set starting k.est to be the maximum number of components
  k.est <- min(dim(phenotype.mx)[1],dim(phenotype.mx)[2])
  cat("Initial k.est = [",k.est,"]\n")

  iter <- 1
  repeat{
    ica.list <- list()

    # 1st run of ICA
    ica.list[[1]] <- fastICA::fastICA(phenotype.mx, k.est,
                             alg.typ = "parallel",method = "C",
                             fun = "logcosh" ,                            # function that should be used to estimate ICs, default is logcosh
                             alpha = 1,row.norm = FALSE,                  # row.norm is set to false since the phenotype.mx is scaled separately
                             maxit=500,tol = 0.0001, verbose = FALSE)
    # 2nd run of ICA
    ica.list[[2]] <- fastICA::fastICA(phenotype.mx, k.est,
                             alg.typ = "parallel",method = "C",
                             fun = "logcosh" ,                            # function that should be used to estimate ICs, default is logcosh
                             alpha = 1,row.norm = FALSE,                  # row.norm is set to false since the phenotype.mx is scaled separately
                             maxit=500,tol = 0.0001, verbose = FALSE)

    # 3rd run of ICA
    #      ica.list[[3]] <- fastICA(t(phenotype.mx), k.est,
    #                               alg.typ = "parallel",method = "C",
    #                               fun = "logcosh" ,                            # function that should be used to estimate ICs, default is logcosh
    #                               alpha = 1,row.norm = FALSE,                  # row.norm is set to false since the phenotype.mx is scaled separately
    #                               maxit=500,tol = 0.0001, verbose = FALSE)
    # calculate the correlation between the signals
    corr.mx1 <- cor(ica.list[[1]]$S, ica.list[[2]]$S)
    #      corr.mx2 <- cor(ica.list[[1]]$S, ica.list[[3]]$S)
    #      corr.mx3 <- cor(ica.list[[2]]$S, ica.list[[3]]$S)

    # The number of pairs of signals having an absolute correlation higher than 0.5 will be used as the next k
    # Meaning = number of replicating signals
    #      k.est.update <- round(mean(sum(abs(corr.mx1) > 0.5), sum(abs(corr.mx2) > 0.5),sum(abs(corr.mx3) > 0.5)))
    k.est.update <- sum(abs(corr.mx1) > 0.5)


    # save the information for the current run
    summary.vec <- rbind(summary.vec, c("run" = run.index,"iter" = iter,
                                        "k.est"= k.est,"rep.ratio"= k.est.update/k.est))

    if(k.est.update/k.est <= 1 & abs((k.est.update /k.est)-1) < 0.01){
      k.vec <- k.est
      cat("Final Value of k in run #", run.index,"= [",k.vec,"] \n")
      break;
    }

    if (iter == max.iter){
      cat("Algorithm did not converge after reaching maximum iterations = ",max.iter,"\n")
      k.vec <- k.est
      cat("Last estimate of k in run # ",run.index," = [",k.vec,"] \n")
      break;

    }

    k.est <- k.est.update
    cat("Run #",run.index,"iteration =",iter,"Updated k.est = [",k.est,"]\n")
    iter <- iter + 1

  }

  return(list("k.est" =k.vec,"summary.vec" = summary.vec ))

}
