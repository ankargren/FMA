PIE <- function(y, X, rhovec = seq(0.2, 0.8, by = 0.2), nreps = 100, avertype = "all", estimators = c("AIC", "BIC", "MMA", "JMA")){
  if (!(is.matrix(y) == TRUE || dim(y)[2] != 1)) {
    warning("y must be a n x 1 matrix.")
    stop
  }
  s.matrix <- Generate.S(Fix.X = 1, Potential.X = ncol(X) - 1, avertype = avertype)
  n <- nrow(X)
  allMods <- EstAllModels(X = X, Xnew = X, y = y, s = s.matrix)
  sigma2 <- allMods$SEs[length(allMods$SEs)]
  originalpreds <- sapply(estimators, function(x) {allMods$preds %*% FMA(x, allMods)})
  temp <- function(y, X, avertype, rho) {
    yperturb <- y + rnorm(nrow(y), 0, rho)
    allModsperturb <- EstAllModels(X = X, Xnew = X, y = yperturb, s = s.matrix)
    perturbpreds <- sapply(estimators, function(x) {allModsperturb$preds %*% FMA(x, allModsperturb)})
    return(sqrt(colMeans((perturbpreds - originalpreds)^2)/sigma2))
  }
  returnval <- sapply(rhovec, function(rho) rowMeans(replicate(nreps, temp(y, X, avertype, rho))))
  colnames(returnval) <- rhovec
  rownames(returnval) <- estimators
  return(returnval)
}