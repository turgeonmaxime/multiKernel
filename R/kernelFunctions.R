# Define kernels----
#' Kernel functions
#' 
#' @param X n x q matrix of covariates
#' @param X_test matrix of test covariates
#' @return Gram matrix for specified kernel
#' @export
#' @importFrom kernlab kernelMatrix
linearKernel <- function(X, X_test = X) {	
  tcrossprod(X_test, X) 
}

### Quadratic kernel
#' @rdname linearKernel
#' @export 
quadraticKernel <- function(X, X_test = X) { 
  (1 + tcrossprod(X_test, X))^2
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
#' @importFrom kernlab rbfdot
gaussKernel <- function(X, X_test = X, rho = 1) {
  # Kmat <- matrix(NA, nrow=nrow(X), ncol=nrow(X))
  # 
  # for (i in 1:nrow(X)){
  #   Kmat[i,] <- diag(tcrossprod(X_test - X_test[i,],
  #                               X - X[i,]))
  #   
  # }
  # 
  # return(exp(-Kmat/rho))
  kernelMatrix(rbfdot(rho), X_test, X)@.Data
}
