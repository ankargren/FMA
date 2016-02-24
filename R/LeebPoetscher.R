LeebPoetscher <- function(samplesize, gam, setup) {
  simData <- list(X = NULL, y = NULL, mu = NULL, thetan = NULL)
  
  if (setup == 1) {
    eta <- c(0, 0, 1, 1, 0, 1, 1, 1)
  }
  if (setup == 2) {
    eta <- c(0, 0, 1, 1, 0, 0, 0, 0)
  }
  if (setup == 3) {
    eta <- c(0, 0, 1, 1, 0, 1/10, 1/10, 1/10)
  }

  theta0 <- c(3, 1.5, 0, 0, 2, 0, 0, 0)
  sigma  <- 1

  thetan <- theta0 + gam / sqrt(samplesize) * eta
  
  XCov <- outer(1:8, 1:8, function(x, y) 0.5^(abs(x-y)))
  simData$X <- matrix(rnorm(samplesize * 8), ncol = 8) %*% t(chol(XCov))
  simData$mu <- simData$X %*% thetan
  simData$y <- simData$mu + rnorm(samplesize, mean = 0, sd = sigma)
  simData$thetan <- thetan
  return(simData)
}