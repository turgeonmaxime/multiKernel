#' Multivariate kernel regression
#' 
#' This function performs multivariate kernel regression by minimizing a specific loss function
#' 
#' @param response matrix of response variables
#' @param covariate matrix of covariate variables, which are included in the kernel.
#' @param confounder matrix or data.frame of confounder variables, which are not included in the kernel.
#' @param kernel Type of kernel to use.
#' @param tau Tuning parameter.
#' @param pure Logical. Use the pure R version?
#' @return Returns the kernel predictor.
fitMultiKernel <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), tau, pure=FALSE) {
  n <- nrow(response); p <- ncol(response)
  if(is.null(confounder)) Z_mat <- matrix(1, nrow=n, ncol=1) else Z_mat <- model.matrix(~., as.data.frame(confounder))
  
  # alpha_mat <- matrix(0, nrow = n, ncol = p)
  # B <- matrix(0, nrow=ncol(Z_mat), ncol=p)
  
  # Kmat <- NULL
  kernel <- match.arg(kernel)
  if(kernel == "linear") Kmat <- linearKernel(covariate)
  if(kernel == "quadratic") Kmat <- quadraticKernel(covariate)
  if(kernel == "gaussian") Kmat <- gaussKernel(covariate)
  # if(is.null(Kmat)) stop("The requested kernel has not been implemented.")
  
  if (pure) {
    weight_mat <- solve(Kmat + diag(tau, ncol = n, nrow = n))
    
    B <- solve(crossprod(Z_mat, weight_mat %*% Z_mat)) %*% crossprod(Z_mat, weight_mat) %*% response
    alpha_mat <- weight_mat %*% (response - Z_mat %*% B)
    
    # currentLS <- 1
    # newLS <- 0
    counter <- 0
    # while(counter < 10000 && abs(currentLS - newLS) > 1e-8) {
    #   counter <- counter + 1
    #   currentLS <- newLS
    #   
    #   # Update B
    #   B <- lm.fit(x = Z_mat, y = (response - K %*% alpha_mat))$coefficients
    #   
    #   # Update alpha_mat
    #   for (j in 1:p) {
    #     alpha_mat[, j] <- weight_mat %*% (response[, j] - Z_mat %*% B[, j])
    #   }
    #   
    #   newLS <- computeLeastSq(response, Kmat, alpha_mat, Z_mat, B)
    #   
    # }
    out <- list(alpha = alpha_mat, B = B, iter = counter)
    
  } else {
    out <- multiKernel_cpp(response, Z_mat, Kmat, tau)
  }
  
  return(out)
}

#' Cross-validation for multivariate kernel regression
#' 
#' This function performs cross-validation for multivariate kernel regression
#' 
#' @param response matrix of response variables
#' @param covariate matrix of covariate variables, which are included in the kernel.
#' @param confounder matrix or data.frame of confounder variables, which are not included in the kernel.
#' @param kernel Type of kernel to use.
#' @param tau Tuning parameter.
#' @param K number of folds for cross-validation.
#' @param pure Logical. Use the pure R version?
#' @return Returns a list of kernel predictors, indexed by the different values of tau.
cvMultiKernel <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), tau, K = 5, pure = FALSE) {
  n <- nrow(response); p <- ncol(response)
  if(is.null(confounder)) Z_mat <- matrix(1, nrow=n, ncol=1) else Z_mat <- model.matrix(~., as.data.frame(confounder))
  
  # K_mat <- NULL
  kernel <- match.arg(kernel)
  if(kernel == "linear") compKernel <- linearKernel
  if(kernel == "quadratic") compKernel <- quadraticKernel
  if(kernel == "gaussian") compKernel <- gaussKernel
  # if(is.null(Kmat)) stop("The requested kernel has not been implemented.")
  
  # create folds for CV
  folds <- createFold(response, K)
  
  predErr_vect <- rep_len(NA, K)
  
  for (i in 1:K) {
    # Training data
    response_train <- response[unlist(folds[-i]),]
    covariate_train <- covariate[unlist(folds[-i]),]
    Z_mat_train <- Z_mat[unlist(folds[-i]),,drop=FALSE]
    
    Kmat_train <- compKernel(covariate_train)
    
    if(pure) {
      n_train <- nrow(response_train)
      weight_mat <- solve(Kmat_train + diag(tau, ncol = n_train, nrow = n_train))
      
      B <- solve(crossprod(Z_mat_train, weight_mat %*% Z_mat_train)) %*% crossprod(Z_mat_train, weight_mat) %*% response_train
      alpha_mat <- weight_mat %*% (response_train - Z_mat_train %*% B)
      out <- list(alpha = alpha_mat, B = B)
    } else {
      out <- multiKernel_cpp(response_train, Z_mat_train, Kmat_train, tau)
    }
      
    # compute prediction error on test data
    response_test <- response[unlist(folds[i]), ]
    covariate_test <- covariate[unlist(folds[i]), ]
    Z_mat_test <- Z_mat[unlist(folds[i]), , drop=FALSE]
    
    Kmat_test <- compKernel(covariate_train, covariate_test)
    prediction <- crossprod(Kmat_test, out$alpha) + Z_mat_test %*% out$B
    
    predErr_vect[i] <- norm(response_test - prediction, type = "F")^2
  }
  
  return(list(predErr = mean(predErr_vect), predErr_vect))
}

#' Selection of tuning parameter for multivariate kernel regression
#' 
#' This function performs cross-validation for multivariate kernel regression
#' 
#' @param response matrix of response variables
#' @param covariate matrix of covariate variables, which are included in the kernel.
#' @param confounder matrix or data.frame of confounder variables, which are not included in the kernel.
#' @param kernel Type of kernel to use.
#' @param tau_seq Sequence of tuning parameters.
#' @return Returns a list of kernel predictors, indexed by the different values of tau.
selectMultiKernel <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), tau_seq) {
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
  
  output_list <- vector("list", length(tau_seq))
  names(output_list) <- tau_seq
  for (tau in tau_seq) {
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
    
    output_list[[which(tau_seq == tau)]] <- list(alpha = alpha_mat, 
                                                 B = B, 
                                                 iter = counter,
                                                 LS = newLS,
                                                 BIC = 2 * newLS + log(n) * (sum(B > 0) + sum(alpha_mat > 0)))
  }
  
  return(output_list)
}

selectMultiKernel2 <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), tau_seq) {
  n <- nrow(response); p <- ncol(response)
  if(is.null(confounder)) Z_mat <- matrix(1, nrow=n, ncol=1) else Z_mat <- model.matrix(~., as.data.frame(confounder))
  
  # K <- NULL
  kernel <- match.arg(kernel)
  if(kernel == "linear") K <- linearKernel(X)
  if(kernel == "quadratic") K <- quadraticKernel(X)
  if(kernel == "gaussian") K <- gaussKernel(X)
  # if(is.null(K)) stop("The requested kernel has not been implemented.")
  
  output_list <- vector("list", length(tau_seq))
  names(output_list) <- tau_seq
  for (tau in tau_seq) {
    out <- multiKernel_cpp(response, Z_mat, K, tau)
    output_list[[which(tau_seq == tau)]] <- out
  }
  
  return(output_list)
}
