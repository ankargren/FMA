Hansen <- function(alpha, n, r2, J) {
  const <- sqrt(r2 / (2 * alpha * zeta(1 + 2 * alpha) - zeta(1 + 2 * alpha) * r2))
  Theta <- const * sqrt(2*alpha)*(1:J)^(-alpha-1/2)
  X <- cbind(1, matrix(rnorm((J-1)*n), nrow = n))
  mu <- X %*% Theta
  y <- as.matrix(mu + rnorm(n))
  return(list(X = X, y = y, mu = mu))
}