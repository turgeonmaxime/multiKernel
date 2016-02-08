# Try out the algorithm on a toy example----
p <- 4; q <- 10; n <- 50
Y <- matrix(rnorm(p * n), nrow = n)
X <- matrix(rnorm(q * n), nrow = n)
Z <- matrix(1, nrow = n, ncol = 1)

Y <- Y + X %*% matrix(1.5, nrow=q, ncol=p) + 2

fitMultiKernel(Y, X, tau = 1)
fitMultiKernel2(Y, X, tau = 1)

compare <- microbenchmark(fitMultiKernel(Y, X, tau = 1),
                          fitMultiKernel2(Y, X, tau = 1),
                          times = 1000)
compare 