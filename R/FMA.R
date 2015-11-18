FMA <- function(X, y, weighting.method, submodels, Xnew, allMods) {
  if (missing(weighting.method) == TRUE) {
    stop("Need to specify the method to use for weighting.")
  }
  
  if (missing(submodels) == TRUE) {
    stop("Need to specify which models to use for weighting.")
  }
  
  if ((missing(X) | missing(y)) == TRUE) {
    stop("Need to supply data.")
  }
  if (missing(Xnew) == TRUE) {
    predictflag <- 1
  } else {
    predictflag <- 0
  }
  
  if (missing(allMods) == TRUE) {
    estimateflag <- 1
  } else {
    estimateflag <- 0
  }
  
  if (!(weighting.method %in% c("AIC", "BIC", "JMA", "MMA"))) {
    stop("The weighting method is incorrectly specified.")
  }
  
  
}