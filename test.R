# Try out the algorithm on a toy example----
p <- 100; q <- 10; n <- 500
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
cvMultiKernel(Y, X, tau = 1, pure=TRUE)
set.seed(12345)
cvMultiKernel(Y, X, tau = 1, pure=FALSE)

# Selection of optimal tuning parameter
foo1 <- selectMultiKernel(Y, X, tau_seq = seq(0.1, 1, length.out = 10), pure=TRUE)
foo2 <- selectMultiKernel(Y, X, tau_seq = seq(0.1, 1, length.out = 10), pure=FALSE)
