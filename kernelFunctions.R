# Define kernels----
#' Kernel functions
#' 
#' @param X n x q matrix of covariates
#' @return Gram matrix for specified kernel
linearKernel <- function(X) {	
  tcrossprod(X) 
}

### Quadratic kernel
#' @rdname linearKernel
quadraticKernel <- function(X) { 
  (1 + tcrossprod(X))^2
}

### Identical by state (IBS) kernel
#' @rdname linearKernel
IBSkernel <- function(X) {
  n <- nrow(X)
  K <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in i:n) {
      K[i,j] <- K[j,i] <- sum(2 * (X[j, ] == X[i, ]) + 
                                (abs(X[j, ] - X[i, ]) == 1))
    }
  }
  
  K <- K/ncol(X)
  return(K)
}

### Gaussian kernel
#' @param rho scaling parameter
#' @rdname linearKernel
gaussKernel <- function(X, rho = 1) {
  Kmat <- matrix(NA, nrow=nrow(X), ncol=nrow(X))
  
  for (i in 1:nrow(X)){
    Kmat[i,] <- diag(tcrossprod(X-X[i,]))
    
  }
  
  return(exp(-Kmat/rho))
}
