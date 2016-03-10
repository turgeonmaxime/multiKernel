predict.multiKernel <- function(object, covariate, confounder=NULL) {
  compKernel <- switch(object$kernel,
                       linear = linearKernel,
                       quadratic = quadraticKernel,
                       gaussian = gaussKernel,
                       stop("The requested kernel has not been implemented.", 
                            call. = FALSE))
  n <- nrow(covariate)
  if(is.null(confounder)) {
    Z_mat <- matrix(1, nrow=n, ncol=1)
    } else Z_mat <- model.matrix(~., as.data.frame(confounder))
  Kmat_test <- compKernel(object$covTrain, covariate)
  if (object$intercept) {
    fitted <- (Kmat_test %*% object$alpha) + Z_mat %*% object$B
  } else {
    fitted <- Kmat_test %*% object$alpha
  }
  return(fitted)
}