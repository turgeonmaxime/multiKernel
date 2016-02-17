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
#' @return Returns the estimated prediction error, as well as a separate value for each fold
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
#' and selects the optimal tuning parameter among a user-specified collection
#' 
#' @param response matrix of response variables
#' @param covariate matrix of covariate variables, which are included in the
#'   kernel.
#' @param confounder matrix or data.frame of confounder variables, which are not
#'   included in the kernel.
#' @param kernel Type of kernel to use.
#' @param tau_seq Sequence of tuning parameters.
#' @param K number of folds for cross-validation.
#' @param pure Logical. Use the pure R version?
#' @return Returns a list of kernel predictors, indexed by the different values
#'   of tau.
selectMultiKernel <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), tau_seq, K = 5, pure = FALSE) {
  
  predError_vect <- vector("numeric", length(tau_seq))
  names(predError_vect) <- tau_seq
  for (tau in tau_seq) {
    out <- cvMultiKernel(response, covariate, confounder, kernel, tau, K, pure)
    predError_vect[as.character(tau)] <- out$predErr
  }
  
  tau_opt <- tau_seq[which.min(predError_vect)]
  out <- fitMultiKernel(response, covariate, confounder, kernel, tau_opt, pure)
  
  return(list(tau_opt, predError_vect, out))
}