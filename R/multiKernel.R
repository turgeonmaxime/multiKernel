#' Multivariate kernel-machine regression
#' 
#' This function performs multivariate kernel-machine regression by minimizing a
#' specific loss function
#' 
#' If \code{confounder = NULL}, \code{intercept = FALSE}, and \code{response}
#' contains only one response variable, then this is equivalent to kernel ridge
#' regression.
#' 
#' @param response matrix of response variables
#' @param covariate matrix of covariate variables, which are included in the
#'   kernel.
#' @param confounder matrix or data.frame of confounder variables, which are not
#'   included in the kernel.
#' @param kernel Type of kernel to use.
#' @param intercept Should we include an intercept?
#' @param tau Tuning parameter.
#' @param pure Logical. Use the pure R version?
#' @param ... Extra parameters to be passed to the kernel function.
#' @return Returns the kernel predictor.
#' @export
#' @useDynLib multiKernel
#' @importFrom Rcpp sourceCpp
fitMultiKernel <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), intercept = TRUE,
                           tau, pure = FALSE, ...) {
  response <- as.matrix(response)
  n <- nrow(response); p <- ncol(response)
  if(is.null(confounder)) Z_mat <- matrix(1, nrow=n, ncol=1) else Z_mat <- model.matrix(~., as.data.frame(confounder))
  
  kernel <- match.arg(kernel)
  Kmat <- switch(kernel,
                 linear = linearKernel(covariate),
                 quadratic = quadraticKernel(covariate),
                 gaussian = gaussKernel(covariate, ...),
                 stop("The requested kernel has not been implemented.", 
                      call. = FALSE))
  
  if (pure) {
    weight_mat <- solve(Kmat + diag(tau, ncol = n, nrow = n))
    
    if (intercept || !is.null(confounder)) {
      B <- solve(crossprod(Z_mat, weight_mat %*% Z_mat)) %*% crossprod(Z_mat, weight_mat) %*% response
      alpha_mat <- weight_mat %*% (response - Z_mat %*% B)
      out <- list(alpha = alpha_mat, B = B)
    } else {
      alpha_mat <- weight_mat %*% response
      out <- list(alpha = alpha_mat)
    }
    
  } else {
    if (intercept || !is.null(confounder)) {
      out <- multiKernel_cpp(response, Z_mat, Kmat, tau)
    } else {
      out <- multiKernel_noCon_cpp(response, Kmat, tau)
    }
  }
  out$kernel <- kernel
  out$intercept <- intercept
  out$covTrain <- covariate
  class(out) <- "multiKernel"
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
#' @param intercept Should we include an intercept?
#' @param tau Tuning parameter.
#' @param K number of folds for cross-validation.
#' @param pure Logical. Use the pure R version?
#' @param ... Extra parameters to be passed to the kernel function.
#' @return Returns the estimated prediction error, as well as a separate value for each fold
#' @export
cvMultiKernel <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), intercept = TRUE, 
                          tau, K = 5, pure = FALSE, ...) {
  n <- nrow(response); p <- ncol(response)
  # if(is.null(confounder)) Z_mat <- matrix(1, nrow=n, ncol=1) else Z_mat <- model.matrix(~., as.data.frame(confounder))
  
  kernel <- match.arg(kernel)
  compKernel <- switch(kernel,
                       linear = linearKernel,
                       quadratic = quadraticKernel,
                       gaussian = function(t) gaussKernel(t, ...),
                       stop("The requested kernel has not been implemented.", 
                            call. = FALSE))
  
  # create folds for CV
  folds <- createFold(response, K)
  fitted_values <- matrix(NA, nrow=n, ncol=p)
  
  for (i in 1:K) {
    # Training data
    response_train <- response[unlist(folds[-i]),,drop=FALSE]
    covariate_train <- covariate[unlist(folds[-i]),]
    confounder_train <- confounder[unlist(folds[-i]),,drop=FALSE]
    
    # Kmat_train <- compKernel(covariate_train)
    out <- fitMultiKernel(response_train, covariate_train, confounder_train, kernel, intercept, tau, pure, ...)
    
    # if(pure) {
    #   n_train <- nrow(response_train)
    #   weight_mat <- solve(Kmat_train + diag(tau, ncol = n_train, nrow = n_train))
    #   
    #   B <- solve(crossprod(Z_mat_train, weight_mat %*% Z_mat_train)) %*% crossprod(Z_mat_train, weight_mat) %*% response_train
    #   alpha_mat <- weight_mat %*% (response_train - Z_mat_train %*% B)
    #   out <- list(alpha = alpha_mat, B = B)
    # } else {
    #   out <- multiKernel_cpp(response_train, Z_mat_train, Kmat_train, tau)
    # }
      
    # compute prediction error on test data
    response_test <- response[unlist(folds[i]), ]
    covariate_test <- covariate[unlist(folds[i]), ]
    confounder_test <- confounder[unlist(folds[i]), , drop=FALSE]
    
    # Kmat_test <- compKernel(covariate_train, covariate_test)
    # fitted_values[unlist(folds[i]),] <- (Kmat_test %*% out$alpha) + Z_mat_test %*% out$B
    fitted_values[unlist(folds[i]),] <- predict(out, covariate_test, confounder_test)
  }
  
  return(list(predErr = mean((response - fitted_values)^2), fitted_values=fitted_values))
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
#' @param intercept Should we include an intercept?
#' @param tau_seq Sequence of tuning parameters.
#' @param K number of folds for cross-validation.
#' @param pure Logical. Use the pure R version?
#' @param ... Extra parameters to be passed to the kernel function.
#' @return Returns a list of kernel predictors, indexed by the different values
#'   of tau.
#' @export
selectMultiKernel <- function(response, covariate, confounder = NULL, kernel = c("linear", "quadratic", "gaussian"), intercept = TRUE, 
                              tau_seq, K = 5, pure = FALSE, ...) {
  
  predError_vect <- vector("numeric", length(tau_seq))
  fitted_list <- vector("list", length(tau_seq))
  names(predError_vect) <- names(fitted_list) <- tau_seq
  
  for (tau in tau_seq) {
    out <- cvMultiKernel(response, covariate, confounder, kernel, intercept, tau, K, pure, ...)
    predError_vect[as.character(tau)] <- out$predErr
    fitted_list[[as.character(tau)]] <- out$fitted_values
  }
  
  tau_opt <- tau_seq[which.min(predError_vect)]
  out <- fitMultiKernel(response, covariate, confounder, kernel, intercept, tau_opt, pure, ...)
  
  return(list(tau_opt, predError_vect, fitted_list, out))
}