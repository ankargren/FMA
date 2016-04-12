Hansen <- function(alpha, n, r2, J, JJ = 1000, rho = 0, equi = FALSE) {
  const <- sqrt(r2 / (2 * alpha * (zeta(1 + 2 * alpha) - zeta(1 + 2 * alpha) * r2)))
  Theta <- const * sqrt(2*alpha)*(1:J)^(-alpha-1/2)
  
  # Compute variance of omitted variables term
  if (rho != 0) {
    thetatemp <- const * sqrt(2 * alpha)/(J:JJ)^(alpha + 0.5)
    lambda <- sum(outer(thetatemp, thetatemp, FUN = function(x, y) x*y)*
                  outer(J:JJ, J:JJ, FUN = function(x, y) rho^(abs(x-y))))
    tau <- const * sqrt(2 * alpha) * (polylog(z = rho, s = alpha + 0.5) - sum(rho^(1:J)/((1:J)^(0.5+alpha))))
    
    if (equi == TRUE) {
      if (alpha <= 0.5) {
        stop("The decay is too slow. Increase alpha.")
      }
      lambda <- 2 * c^2 * alpha * rho^2 * (zeta(1 + 2*alpha) - sum(1/(1:J)^(alpha + 0.5)))^2
      tau <- sqrt(lambda)
      Sigma <- rbind(matrix(rho, J - 1, J - 1) + diag(1 - rho, J - 1), rep(tau, J - 1),
                     c(rep(tau, J - 1), lambda))
    } else {
      # if not equi
      Sigma <- rbind(cbind(outer(1:(J-1), 1:(J-1), function(x, y) rho^(abs(x-y))), tau * 1/(rho^(1:(J-1)))),
                     c(tau * 1/(rho^(1:(J-1))), lambda))
    }
    
  }

  Z <- cbind(1, rmvn(n, rep(0, J), Sigma))
  X <- Z[, -ncol(Z)]
  mu <- Z %*% c(Theta, 1)
  y <- as.matrix(mu + rnorm(n))
  return(list(X = X, y = y, mu = mu))
}