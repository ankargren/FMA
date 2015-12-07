Generate.S <- function(Fix.X, Potential.X, avertype = "all") {
  if (!(avertype %in% c("all", "nested", "singleton"))) {
    stop("Need to specify the type of submodels to average over.")
  }
  
  n.variables <- Fix.X + Potential.X
  if (avertype == "all") {
    n.submodels <- 2^Potential.X
  } else {
    n.submodels <- Potential.X + ifelse(Fix.X > 0, 1, 0)
  }
  
  s <- matrix(0, nrow = n.submodels, ncol = n.variables)
  if (Fix.X > 0) {
    s[, 1:Fix.X] <- 1
  } 
  
  if (avertype == "all") {
    s[, Potential.X] <- as.matrix(expand.grid(rep(list(c(0,1)), Potential.X)))
  } else if (avertype == "nested") {
    s[ifelse(Fix.X > 0, 2, 1):nrow(s), (1 + Fix.X):n.variables] <- lower.tri(diag(Potential.X), 1)*1
  } else {
    s[ifelse(Fix.X > 0, 2, 1):nrow(s), (1 + Fix.X):n.variables] <- diag(Potential.X)
  }
  s
}
