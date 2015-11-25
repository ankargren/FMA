FMA <- function(weighting.method, allMods, X, y, submodels, Xnew, include.intercept = 1) {
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
    } else {
      X <- as.matrix(X)
      y <- as.matrix(y)
    }
    
    if (submodels %in% c("all", "nested", "singleton")) {
      Fix.Z <- NULL
      Potential.Z <- NULL
      if ((include.intercept == 1) || all(X[, 1] == 1)) {
        Fix.X <- 1
        Potential.X <- 2:ncol(X)
      } else {
        Fix.X <- NULL
        Potential.X <- 1:ncol(X)
      }
      s <- Generate.S(Fix.X, Fix.Z, Potential.X, Potential.Z, avertype = "submodels")
    }
    
    if (is.matrix(submodels) && !(dim(submodels)[2] == dim(X)[2])) {
      stop("Dimensions of the submodel matrix are wrong.")
    } else {
      s <- submodels
    }
  } else {
    estimateflag <- 0
  }
  if (missing(Xnew) == TRUE) {
    predictflag <- 1
  } else {
    predictflag <- 0
  }
  
  if (!(weighting.method %in% c("AIC", "BIC", "JMA", "MMA"))) {
    stop("The weighting method is incorrectly specified.")
  }
  
  if (estimateflag == 1) {
    if (predictflag == 0) {
      Xnew <- matrix(0, 1, ncol(allMods$s))
    }
    allMods <- EstAllModels(X, Xnew, y, s)
  }
  
  if (weighting.method %in% c("MMA", "JMA")) {
    M <- nrow(s)
    if (weighting.method == "MMA") {
      cc <- allMods$SEs[M] * allMods$K
      HH <- crossprod(allMods$etilde)
    } else {
      cc <- matrix(0, nrow = M, ncol = 1)
      HH <- crossprod(allMods$eJMA)
    }
    
    QP <- ipop(c = cc, H = HH, A = rep(1, M), b = 1,
                   l = rep(0, M), u = rep(1, M), r = 0, maxiter = 100)
  }
  
  if (weighting.method %in% c("AIC", "BIC")) {
    if (weighting.method == "AIC") {
      ICvec <- allModels$AIC
    } else {
      ICvec <- allModels$BIC
    }
    maxmod <- max(ICvec)
    ICw <- exp(-1/2*(ICvec-maxmod))/sum(exp(-1/2*(ICvec-maxmod)))
  }
}