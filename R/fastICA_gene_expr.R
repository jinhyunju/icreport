#' ICA analysis based on the package fastICA
#'
#' Minor tweak to fastICA to make it faster
#'
#' @param phenotype.mx Phenotype matrix with diemnsions g x N   \code{phenotype.mx}
#' @param k.est Number of components to be estimated or method to estimate it. \code{k.est}
#' @param scale.pheno Logical value specifying the scaling of row of the phenotype.mx. \code{scale.pheno}
#' @return List with the following entries.
#' @keywords keywords
#'
#' @import MASS
#' @export
#'
#' @examples
#' R code here showing how your function works
fastICA_gene_expr <-function(X, n.comp, fun = "logcosh", alpha = 1,
                             scale.pheno = FALSE, maxit = 200, tol = 1e-04, verbose = TRUE, w.init=NULL) {
    dd <- dim(X)       # dimensions g x N
    d <- dd[dd != 1L]
    if (length(d) != 2L)
        stop("data must be matrix-conformal")
    X <- if (length(d) != length(dd)) matrix(X, d[1L], d[2L])
    else as.matrix(X)

    if (alpha < 1 || alpha > 2)
        stop("alpha must be in range [1,2]")
#    method <- match.arg(method)
#    alg.typ <- match.arg(alg.typ)
#    fun <- match.arg(fun)
    n <- nrow(X)         # g
    p <- ncol(X)         # N

    if (n.comp > min(n, p)) {
        message("'n.comp' is too large: reset to ", min(n, p))
        n.comp <- min(n, p)
    }
    if(is.null(w.init))
        w.init <- matrix(rnorm(n.comp^2),n.comp,n.comp)
    else {
        if(!is.matrix(w.init) || length(w.init) != (n.comp^2))
            stop("w.init is not a matrix or is the wrong size")
    }

    pca.X <- prcomp(X)         # X is still g x N

    K.comp <- pca.X$rotation[,c(1:n.comp)]        # k principal components
    Diag <- diag(c(1/pca.X$sdev[c(1:n.comp)]))    # k standard deviation diagonal matrix

    K <- K.comp %*% Diag

    X1 <- X %*% K   # g x N %*% N x k %*% k %*% k

    X1 <- t(X1)
    X <- t(X)
    ### part where ICA is done
    a <- ica_R_par(X1, n.comp, tol = tol, fun = fun, alpha = alpha, maxit = maxit, verbose = verbose, w.init = w.init)

    # reconstructing data before output
    w <- a %*% t(K)    # k x k %*% k x N = k x N
    S <- w %*% X    # k x N %*% N x g = k x g
    A <- t(w) %*% ginv(w %*% t(w)) # N x k x k x k  = N x k

    X <- t(X) # g X N
    W <- t(a) # N x k
    A <- t(A) # k x N
    S <- t(S) # g x k
    return(list(X = X, K = t(K), W = a, A = A, S = S))
}
