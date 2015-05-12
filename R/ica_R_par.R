#' ICA parameter estimation function
#'
#' Adopted from fastICA package
#'
#' @param X Phenotype matrix with diemnsions g x N
#' @param n.comp Number of components to be estimated or method to estimate it.
#' @param fun Function to be used in ICA estimation
#' @param alpha alpha for ICA, should be in range [1,2].
#' @param maxit Maximum iterations
#' @param tol Threshold for convergence
#' @param verbose If TRUE details of the estimation process are shown.
#' @param w.init Initial value for W, if left unspecified random numbers will be used.
#' @return List with the following entries.
#' @keywords keywords
#'
#' @export
#'
ica_R_par <- function (X, n.comp, tol, fun, alpha, maxit, verbose, w.init) {
  Diag <- function(d) if(length(d) > 1L) diag(d) else as.matrix(d)
  n <- nrow(X)
  p <- ncol(X)
  W <- w.init
  sW <- La.svd(W)
  W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
  W1 <- W
  lim <- rep(1000, maxit)
  it <- 1
  if (fun == "logcosh") {
    if (verbose)
      message("Running FastICA ( logcosh approx. ) \n")
    while (lim[it] > tol && it < maxit) {
      wx <- W %*% X
      gwx <- tanh(alpha * wx)
      v1 <- gwx %*% t(X)/p
      g.wx <- alpha * (1 - (gwx)^2)
      v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      W1 <- v1 - v2
      sW1 <- La.svd(W1)
      W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      W <- W1
      if (verbose)
        message("\r Iteration ", it, " tol = ", format(lim[it + 1]))
      it <- it + 1
    }
  }
  if (fun == "exp") {
    if (verbose)
      message("Symmetric FastICA using exponential approx. to neg-entropy function \n")
    while (lim[it] > tol && it < maxit) {
      wx <- W %*% X
      gwx <- wx * exp(-(wx^2)/2)
      v1 <- gwx %*% t(X)/p
      g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
      v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      W1 <- v1 - v2
      sW1 <- La.svd(W1)
      W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      W <- W1
      if (verbose){
        message("Iteration ", it, " tol = ", format(lim[it + 1]),"\r")
      }

      it <- it + 1
    }
  }
  W
}
