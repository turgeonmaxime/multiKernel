#' Multivariate kernel regression
#' 
#' This function performs multivariate kernel regression by minimizing a specific loss function
#' 
#' @param response matrix of response variables
#' @param covariate matrix of covariate variables, which are included in the kernel.
#' @param confounder matrix or data.frame of confounder variables, which are not included in the kernel.
#' @param kernel Type of kernel to use.
#' @param tau Tuning parameter
#' @return Returns the kernel predictor.
fitMultiKernel <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), tau) {
  n <- nrow(response); p <- ncol(response)
  if(is.null(confounder)) Z_mat <- matrix(1, nrow=n, ncol=1) else Z_mat <- model.matrix(~., as.data.frame(confounder))
  
  alpha_mat <- matrix(0, nrow = n, ncol = p)
  B <- matrix(0, nrow=ncol(Z_mat), ncol=p)
  
  # K <- NULL
  kernel <- match.arg(kernel)
  if(kernel == "linear") K <- linearKernel(X)
  if(kernel == "quadratic") K <- quadraticKernel(X)
  if(kernel == "gaussian") K <- gaussKernel(X)
  # if(is.null(K)) stop("The requested kernel has not been implemented.")
  
  weight_mat <- solve(K + diag(tau, ncol = n, nrow = n))
  
  currentLS <- 1
  newLS <- 0
  counter <- 0
  while(counter < 10000 && abs(currentLS - newLS) > 1e-8) {
    counter <- counter + 1
    currentLS <- newLS
    
    # Update B
    B <- lm.fit(x = Z_mat, y = (response - K %*% alpha_mat))$coefficients
    
    # Update alpha_mat
    for (j in 1:p) {
      alpha_mat[, j] <- weight_mat %*% (response[, j] - Z_mat %*% B[, j])
    }
    
    newLS <- computeLeastSq(response, K, alpha_mat, Z_mat, B)
    
  }
  
  return(list(alpha = alpha_mat, B = B, iter = counter))
}

fitMultiKernel2 <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), tau) {
  n <- nrow(response); p <- ncol(response)
  if(is.null(confounder)) Z_mat <- matrix(1, nrow=n, ncol=1) else Z_mat <- model.matrix(~., as.data.frame(confounder))
  
  # K <- NULL
  kernel <- match.arg(kernel)
  if(kernel == "linear") K <- linearKernel(X)
  if(kernel == "quadratic") K <- quadraticKernel(X)
  if(kernel == "gaussian") K <- gaussKernel(X)
  # if(is.null(K)) stop("The requested kernel has not been implemented.")
  
  out <- multiKernel_cpp(response, Z_mat, K, tau)
  
  return(out)
}

computeLeastSq <- function(Y, K, alpha, Z, B) {
  - 0.5 * norm(Y - K %*% alpha - Z %*% B, type = "F")
}