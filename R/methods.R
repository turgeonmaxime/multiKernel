#' Prediction method for multivariate kernel-machine regression
#' 
#' This method performs prediction for new data, based on the output of a fitted
#' multivariate kernel-machine regression
#' 
#' @param object Object of class \code{multiKernel}, such as the output of the 
#'   function \code{\link{fitMultiKernel}}.
#' @param covariate Matrix of covariate variables, which are included in the 
#'   kernel.
#' @param confounder Matrix or data.frame of confounder variables, which are not
#'   included in the kernel.
#' @param ... Extra parameters to be passed to the kernel function.
#' @export
predict.multiKernel <- function(object, covariate, confounder=NULL, ...) {
  compKernel <- switch(object$kernel,
                       linear = linearKernel,
                       quadratic = quadraticKernel,
                       gaussian = function(t) gaussKernel(t, ...),
                       stop("The requested kernel has not been implemented.", 
                            call. = FALSE))
  n <- nrow(covariate)
  if(is.null(confounder)) {
    Z_mat <- matrix(1, nrow=n, ncol=1)
    } else Z_mat <- model.matrix(~., as.data.frame(confounder))
  Kmat_test <- compKernel(object$covTrain, covariate)
  if (object$intercept || !is.null(confounder)) {
    fitted <- (Kmat_test %*% object$alpha) + Z_mat %*% object$B
  } else {
    fitted <- Kmat_test %*% object$alpha
  }
  return(fitted)
}