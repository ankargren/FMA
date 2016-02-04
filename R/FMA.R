FMA2 <- function (weighting.method, allMods, X, y, submodels, include.intercept = 1, solver = "solve.QP") 
{
  if (missing(weighting.method) == TRUE) {
    stop("Need to specify the method to use for weighting.")
  }
  if (missing(allMods) == TRUE) {
    estimateflag <- 1
    if (missing(submodels) == TRUE) {
      stop("Need to specify which models to use for weighting.")
    }
    if ((missing(X) | missing(y)) == TRUE) {
      stop("Need to supply data.")
    }
    else {
      X <- as.matrix(X)
      y <- as.matrix(y)
    }
    if (is.character(submodels) && submodels %in% c("all", 
                                                    "nested", "singleton")) {
      if ((include.intercept == 1) || all(X[, 1] == 1)) {
        Fix.X <- 1
        Potential.X <- 2:ncol(X)
      }
      else {
        Fix.X <- NULL
        Potential.X <- 1:ncol(X)
      }
      s <- Generate.S(Fix.X, Potential.X, avertype = submodels)
    }
    if (is.matrix(submodels) == TRUE) {
      if (!(dim(submodels)[2] == dim(X)[2])) {
        stop("Dimensions of the submodel matrix are wrong.")
      }
      else {
        s <- submodels
        M <- nrow(s)
      }
    }
  }
  else {
    estimateflag <- 0
    M <- length(allMods$AIC)
  }
  if (!(weighting.method %in% c("AIC", "BIC", "JMA", "MMA"))) {
    stop("The weighting method is incorrectly specified.")
  }
  if (estimateflag == 1) {
    Xnew <- matrix(0, 1, ncol(X))
    allMods <- EstAllModels(X, Xnew, y, s)
  }
  if (weighting.method %in% c("MMA", "JMA")) {
    if (weighting.method == "MMA") {
      cc <- allMods$SEs[M] * allMods$K
      HH <- t(allMods$etilde)
    }
    else {
      cc <- matrix(0, nrow = M, ncol = 1)
      HH <- t(allMods$eJMA)
    }
    if (solver == "LowRankQP") {
      QP <- LowRankQP(Vmat = HH, dvec = cc, Amat = matrix(1, 1, M), bvec = 1, uvec = rep(1, M))
      ICw <- QP$alpha
    }
    if (solver == "solve.QP") {
      a1 <- HH %*% t(HH) + diag(M) * 1e-8
      a3 <- t(rbind(matrix(1, nrow = 1, ncol = M), diag(M), -diag(M)))
      a4 <- rbind(1, matrix(0, nrow = M, ncol = 1), matrix(-1, nrow = M, ncol = 1))
      
      QP <- solve.QP(Dmat = a1, dvec = -as.vector(cc), Amat = a3, bvec = a4, meq = 1)
      ICw <- QP$solution
    }
    if (solver == "ipop") {
      QP <- ipop(c = cc, H = HH %*% t(HH), A = matrix(1, nrow = 1, ncol = M), b = 1, l = matrix(0, nrow = M, ncol = 1),
                 u = matrix(1, nrow = M, ncol = 1), r = 0)
      ICw <- primal(QP)
    }
    
  }
  if (weighting.method %in% c("AIC", "BIC")) {
    if (weighting.method == "AIC") {
      ICvec <- allMods$AIC
    }
    else {
      ICvec <- allMods$BIC
    }
    maxmod <- max(ICvec)
    ICw <- exp(-1/2 * (ICvec - maxmod))/sum(exp(-1/2 * (ICvec - 
                                                          maxmod)))
  }
  as.numeric(ICw)
}