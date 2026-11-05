sim_data <- function(n, bias = 0, hte = 0.1) {
  X1 <- rnorm(n); X3 <- rnorm(n, 1); U <- rnorm(n)
  A <- if (bias == 0) rbinom(n, 1, 0.5) else rbinom(n, 1, plogis(X1 + X3 - 1 + bias * U))
  data.frame(X1 = X1, X3 = X3, A = A,
             Y = 0.5 + 0.2 * A + X3 + bias * U + hte * A * X3 + rnorm(n))
}

sim_bin <- function(n, bias = 0) {
  X <- rnorm(n); U <- rnorm(n)
  A <- if (bias == 0) rbinom(n, 1, 0.5) else rbinom(n, 1, plogis(X + bias * U))
  data.frame(X = X, A = A, Y = rbinom(n, 1, plogis(-0.3 + 0.8 * A + 0.5 * X + 0.3 * A * X + bias * U)))
}
