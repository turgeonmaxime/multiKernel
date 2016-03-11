# Try out the algorithm on a toy example----
p <- 2; q <- 100; n <- 500
Y <- matrix(rnorm(p * n), nrow = n)
X <- matrix(rnorm(q * n), nrow = n)
X_test <- matrix(rnorm(q * n), nrow = n)
Z <- matrix(1, nrow = n, ncol = 1)

Y <- Y + (sin(X)/X) %*% matrix(1.5, nrow=q, ncol=p) + 2
Y_test <- matrix(rnorm(p * n), nrow = n) + (sin(X_test)/X_test) %*% matrix(1.5, nrow=q, ncol=p) + 2

# Rcpp::sourceCpp('multiKernel.cpp')
# source('multiKernel.R')

fit1 <- fitMultiKernel(Y, X, tau = 1, intercept = FALSE)
fit2 <- fitMultiKernel(Y, X, tau = 1, intercept = TRUE)

library(microbenchmark)
compare <- microbenchmark(fitMultiKernel(Y, X, tau = 1),
                          fitMultiKernel2(Y, X, tau = 1),
                          times = 100)
compare

# Do we need while loop?
tau <- 1; K <- linearKernel(X)
beta <- solve(t(Z) %*% solve(K + diag(tau, ncol = n, nrow = n)) %*% Z) %*% t(Z) %*% solve(K + diag(tau, ncol = n, nrow = n)) %*% Y
alpha <- solve(K + diag(tau, ncol = n, nrow = n)) %*% (Y - Z %*% beta)
# No, we don't...

# Cross-validation---
set.seed(12345)
foo1 <- cvMultiKernel(Y, X, tau = 1, pure=FALSE)
foo2 <- cvMultiKernel(Y, X, tau = 1, intercept = FALSE)
set.seed(12345)
cvMultiKernel(Y, X, tau = 1, pure=FALSE)

# Selection of optimal tuning parameter
set.seed(12)
foo1 <- selectMultiKernel(Y, X, tau_seq = seq(0.1, 1, length.out = 10), pure=TRUE)
set.seed(12)
foo2 <- selectMultiKernel(Y, X, tau_seq = seq(0.1, 1, length.out = 10), pure=FALSE)

tau_seq <- seq(1, 10, length.out = 25)
foo <- selectMultiKernel(Y, X, tau_seq = tau_seq, 
                         K=10, pure=FALSE, intercept = FALSE, kernel = "linear")
predErr <- rep_len(NA, length(tau_seq))
for (i in 1:length(tau_seq)) {
  out <- fitMultiKernel(Y, X, intercept = TRUE, kernel = "linear", tau = tau_seq[i], pure=FALSE)
  pred <- predict(out, X_test)
  predErr[i] <-  mean((pred - Y_test)^2)
}

plot(tau_seq, foo[[2]], type='b', pch=19, cex=0.5,
     xlab="tau", ylab="Prediction error")
plot(tau_seq, predErr, type='b', pch=19, cex=0.5,
     xlab="tau", ylab="Prediction error")
lines(lowess(data.frame(tau_seq, foo[[2]])), lwd=2, col='blue')
lines(tau_seq, foo[[2]], type='b', pch=19, cex=0.5, col='blue', lwd=2)
lines(tau_seq, apply(foo[[3]], 2, median), type='b', pch=19, cex=0.5, col='red', lwd=2)

pred_errors <- lapply(foo1[[4]], function(mat) mat - Y)
sapply(pred_errors, function(mat) mean(mat^2))


########################
library(CVST)

ns = noisySinc(100)
nsTest = noisySinc(1000)
# Kernel ridge regression
krr = constructKRRLearner()
m <- m2 <- vector("numeric", 100)
lambda_vect <- seq(0.001, 1, length.out = 100)

index <- 0
n <- getN(ns)
Kmat <- linearKernel(ns$x)
for (lam in lambda_vect) {
  index <- index + 1
  foo <- krr$learn(ns, list(kernel="rbfdot", lambda=lam))
  pred <- krr$predict(foo, nsTest)
  m[index] <- sum((pred - nsTest$y)^2) / getN(nsTest)
  
  out <- fitMultiKernel(as.matrix(ns$y), ns$x, intercept = FALSE, kernel = "gaussian", tau = lam * n, pure=FALSE)
  pred <- predict(out, nsTest$x)
  # pred <- t(foo$alpha) %*% t(K)
  m2[index] <- sum((pred - nsTest$y)^2) / getN(nsTest)
}
plot(m, m2)
plot(x=lambda_vect, y=m, type='b', pch=19, cex=0.5, ylim=range(c(m, m2[-1])))
lines(x=lambda_vect, y=m2, type='b', pch=19, cex=0.5, col='blue')
cv_out <- CV(ns, krr, constructParams(kernel="rbfdot", sigma=100, lambda=lambda_vect), fold = 10, verbose = TRUE)

fit <- fitMultiKernel(as.matrix(ns$y), ns$x, tau = getN(ns) *lambda_vect[1], intercept = FALSE)
K <- linearKernel(ns$x, nsTest$x)
pred <- K %*% fit$alpha
sum((pred - nsTest$y)^2) / getN(nsTest)

foo <- krr$learn(ns, list(kernel="vanilladot", lambda=lambda_vect[1]))
pred <- krr$predict(foo, nsTest)
sum((pred - nsTest$y)^2) / getN(nsTest)