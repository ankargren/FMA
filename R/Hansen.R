Hansen <- function(alpha, n, r2, J, JJ = 1000, rho = 0, equi = FALSE, include.intercept = TRUE) {
  const <- sqrt(r2 / (2 * alpha * (zeta(1 + 2 * alpha) - zeta(1 + 2 * alpha) * r2)))
  Theta <- const * sqrt(2*alpha)*(1:J)^(-alpha-1/2)
  
  if (include.intercept == FALSE) {
    J_random <- J
  } else {
    J_random <- J - 1
  }
  
  # Compute variance of omitted variables term
  if (rho != 0) {
    thetatemp <- const * sqrt(2 * alpha)/(J:JJ)^(alpha + 0.5)
    lambda <- sum(outer(thetatemp, thetatemp, FUN = function(x, y) x*y)*
                  outer(J:JJ, J:JJ, FUN = function(x, y) rho^(abs(x-y))))
    tau <- const * sqrt(2 * alpha) * (polylog(z = rho, s = alpha + 0.5, n.sum = JJ) - sum(rho^(1:J)/((1:J)^(0.5+alpha))))
    
    if (equi == TRUE) {
      if (alpha <= 0.5) {
        stop("The decay is too slow. Increase alpha.")
      }
      lambda <- 2 * c^2 * alpha * rho^2 * (zeta(1 + 2*alpha) - sum(1/(1:J)^(alpha + 0.5)))^2
      tau <- sqrt(lambda)
      Sigma <- rbind(matrix(rho, J_random, J_random) + diag(1 - rho, J_random), rep(tau, J_random),
                     c(rep(tau, J_random), lambda))
    } else {
      # if not equi
      Sigma <- rbind(cbind(outer(1:(J_random), 1:(J_random), function(x, y) rho^(abs(x-y))), tau * 1/(rho^(1:(J_random)))),
                     c(tau * 1/(rho^(1:(J_random))), lambda))
    }
    Z <- rmvn(n, rep(0, J_random + 1), Sigma)
  } else {
    varw <- 2*const^2*alpha*zeta(1+2*alpha) - sum(Theta^2)
    Z <- cbind(matrix(rnorm(n*(J_random)), n, J_random), rnorm(n, 0, sqrt(varw)))
  }
  # Z = (X, w), see the simulation details
  if (include.intercept == TRUE) {
    Z <- cbind(1, Z)
  }
  X <- Z[, -ncol(Z)]
  mu <- Z %*% c(Theta, 1)
  y <- as.matrix(mu + rnorm(n))
  return(list(X = X, y = y, mu = mu))
}