# Try out the algorithm on a toy example----
p <- 100; q <- 10; n <- 500
Y <- matrix(rnorm(p * n), nrow = n)
X <- matrix(rnorm(q * n), nrow = n)
Z <- matrix(1, nrow = n, ncol = 1)

Y <- Y + X %*% matrix(1.5, nrow=q, ncol=p) + 2

Rcpp::sourceCpp('multiKernel.cpp')
source('multiKernel.R')

fitMultiKernel(Y, X, tau = 1)
fitMultiKernel2(Y, X, tau = 1)

compare <- microbenchmark(fitMultiKernel(Y, X, tau = 1),
                          fitMultiKernel2(Y, X, tau = 1),
                          times = 100)
compare 