# Try out the algorithm on a toy example----
p <- 10; q <- 100; n <- 500
Y <- matrix(rnorm(p * n), nrow = n)
X <- matrix(rnorm(q * n), nrow = n)
Z <- matrix(1, nrow = n, ncol = 1)

Y <- Y + X %*% matrix(1.5, nrow=q, ncol=p) + 2

Rcpp::sourceCpp('multiKernel.cpp')
source('multiKernel.R')

fit1 <- fitMultiKernel(Y, X, tau = 1)
fit2 <- fitMultiKernel2(Y, X, tau = 1)

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
foo1 <- cvMultiKernel(Y, X, tau = 1, pure=TRUE)
set.seed(12345)
cvMultiKernel(Y, X, tau = 1, pure=FALSE)

# Selection of optimal tuning parameter
set.seed(12)
foo1 <- selectMultiKernel(Y, X, tau_seq = seq(0.1, 1, length.out = 10), pure=TRUE)
set.seed(12)
foo2 <- selectMultiKernel(Y, X, tau_seq = seq(0.1, 1, length.out = 10), pure=FALSE)

foo <- selectMultiKernel(Y, X, tau_seq = seq(1, 50, length.out = 50), K=10, pure=FALSE)
plot(seq(1, 50, length.out = 50), foo[[2]], type='b', pch=19, cex=0.5,
     xlab="tau", ylab="Prediction error")
boxplot(foo[[3]])
lines(seq(1, 50, length.out = 50), foo[[2]], type='b', pch=19, cex=0.5, col='blue', lwd=2)
lines(seq(1, 50, length.out = 50), apply(foo[[3]], 2, median), type='b', pch=19, cex=0.5, col='red', lwd=2)

pred_errors <- lapply(foo1[[4]], function(mat) mat - Y)
sapply(pred_errors, function(mat) mean(mat^2))
