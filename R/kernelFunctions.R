# Define kernels----
#' Kernel functions
#' 
#' @param X n x q matrix of covariates
#' @param X_test matrix of test covariates
#' @return Gram matrix for specified kernel
#' @export
linearKernel <- function(X, X_test = X) {	
  tcrossprod(X, X_test) 
}

### Quadratic kernel
#' @rdname linearKernel
#' @export 
quadraticKernel <- function(X, X_test = X) { 
  (1 + tcrossprod(X, X_test))^2
}

### Identical by state (IBS) kernel
#' @rdname linearKernel
#' @export 
IBSkernel <- function(X, X_test = X) {
  n <- nrow(X)
  K <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in i:n) {
      K[i,j] <- K[j,i] <- sum(2 * (X[j, ] == X_test[i, ]) + 
                                (abs(X[j, ] - X_test[i, ]) == 1))
    }
  }
  
  K <- K/ncol(X)
  return(K)
}

### Gaussian kernel
#' @param rho scaling parameter
#' @rdname linearKernel
#' @export
gaussKernel <- function(X, X_test = X, rho = 1) {
  Kmat <- matrix(NA, nrow=nrow(X), ncol=nrow(X))
  
  for (i in 1:nrow(X)){
    Kmat[i,] <- diag(tcrossprod(X - X[i,],
                                X_test - X_test[i,]))
    
  }
  
  return(exp(-Kmat/rho))
}
