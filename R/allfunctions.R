Generate.S <- cmpfun(function(Fix.X, Fix.Z, Potential.X, Potential.Z, avertype = "all") {
  if (!(avertype %in% c("all", "nested", "singleton"))) {
    stop("Need to specify the type of submodels to average over.")
  }

  n.variables <- length(c(Fix.X, Fix.Z, Potential.X, Potential.Z))
  if (avertype == "all") {
    n.submodels <- 2^length(c(Potential.X, Potential.Z))
  } else {
    n.submodels <- length(c(Potential.X, Potential.Z))
  }
  n.potential <- length(c(Potential.X, Potential.Z))
  n.X <- length(c(Fix.X, Potential.X))
  n.Z <- length(c(Fix.Z, Potential.Z))
  s <- matrix(NaN, nrow = n.submodels, ncol = n.variables)
  s[, c(Fix.X, n.X+Fix.Z)] <- 1
  if (avertype == "all") {
    s[, c(Potential.X, n.X+Potential.Z)] <- as.matrix(expand.grid(rep(list(c(0,1)), n.potential)))
  } else if (avertype == "nested") {
    s[, c(Potential.X, n.X+Potential.Z)] <- lower.tri(diag(n.potential), 1)*1
  } else {
    s[, c(Potential.X, n.X+Potential.Z)] <- diag(n.potential)
  }
  s
})

Population.Parameters <- function(J.X, J.Z, Cov.Const, Var.Const, Include.Intercept) {
  Covariance.X  <- outer(1:J.X, 1:J.X, FUN = function(x, y) Cov.Const^abs(x-y))
  Covariance.Z  <- outer(1:J.Z, 1:J.Z, FUN = function(x, y) Cov.Const^abs(x-y))
  Covariance.XZ <- outer(1:J.X, 1:J.Z, FUN = function(x, y) Cov.Const^(1 + abs(x-y)))
  
  # Construct a large covariance matrix and scale it so that variances are not 1
  Covariance <- rbind(cbind(Covariance.X, Covariance.XZ), cbind(t(Covariance.XZ), Covariance.Z))
  Covariance <- diag(sqrt(Var.Const), J.X + J.Z) %*% Covariance %*% diag(sqrt(Var.Const), J.X + J.Z) ### The variances of all variables are 2.
  
  # Population mean
  Mu <- ifelse(Include.Intercept,2,0) # Intercept. If Include.Intercept <- FALSE, it is ignored.
  
  # The current simulation setting is Hansen (2007)
  Alpha <- 0.5; Const <- 2 # Alpha and c values in Hansen (2007)
  Theta.X <- rep(NA, J.X); Theta.Z <- rep(NA, J.Z)
  for(j in 1:J.X){
    Theta.X[j] <- Const * sqrt(2 * Alpha) * j^(-Alpha - 1/2)
  }
  for(j in 1:J.Z){
    Theta.Z[j] <- Const * sqrt(2 * Alpha) * j^(-Alpha - 1/2)
  }
  
  Col.Index <- c(Fix.X, J.X + Fix.Z, Potential.X, J.X + Potential.Z)
  
  if (Include.Intercept == 1) {
    s.matrix <- cbind(1, Generate.S(Fix.X, Fix.Z, Potential.X, Potential.Z))
  } else {
    s.matrix <- Generate.S(Fix.X, Fix.Z, Potential.X, Potential.Z)
  }
  
  
  # Zeta.Func <- function(x) sum(sapply( c(1:1e04), function(y) 1/(y^x) ))
  # Pop.R.Square <- 2*Alpha*Zeta.Func(1+2*Alpha)*Const^2/( 1+Zeta.Func(1+2*Alpha)*Const^2 ) # Maybe not right
  return(list(Mu = Mu, Covariance = Covariance, Theta.X = Theta.X, Theta.Z = Theta.Z, Col.Index = Col.Index, S.matrix = s.matrix))
}

Check.Correct.Model <- function(input.model, target.model){
  all(input.model == target.model)
}

Create.FMA <- function(w, name, allMods) {
  assign(name, list(w = as.numeric(w), coefs = as.numeric(allMods$coefficients %*% w),
                    preds = as.numeric(allMods$preds %*% w)), envir = .GlobalEnv)
}

Wrapper.Glmnet <- function(alphaval, penalty.weights, methodname, include.intercept = Include.Intercept,
                           ry = Ry, rxz = Rxz, col.index = Col.Index, new.OBS = New.OBS, new.Ry = New.Ry) {
  cv <- cv.glmnet(x = rxz[, col.index], y = ry, nfolds = 10, alpha = alphaval, nlambda = 100, type.measure = "mse",
                  standardize = TRUE, intercept = include.intercept, penalty.factor = penalty.weights)
  Pred <- predict(cv, newx = as.matrix(new.OBS), s = "lambda.min")
  return(c(mean(Pred - new.Ry), mean((Pred - new.Ry)^2)))
}

Data.Generation <- function(samplesize, covariance, error.sd, mu, theta.x, col.index, j.x, j.z, new.samplesize) {
  # Generate X and Z
  Rxz <- rmvnorm(samplesize, mean = rep(0, nrow(covariance)), sigma = covariance)
  
  # Generate errors and y
  Error <- matrix(rnorm(samplesize, 0, error.sd), ncol = 1)
  Ry <- matrix(mu, nrow = samplesize, ncol = 1) + Rxz %*% matrix(c(theta.x, rep(0, j.z)), ncol = 1) + Error
  
  # We only observe a subset of the predictors
  OBS <- data.frame(cbind(Ry, Rxz[, col.index]))
  colnames(OBS) <- c('y', c(paste('x', 1:j.x, sep=''), paste('z', 1:j.z, sep=''))[col.index])
  rm(Error)
  
  # Generate new observations (the evaluation set)
  New.Rxz <- rmvnorm(new.samplesize, mean = rep(0, nrow(covariance)), sigma = covariance)
  Error <- matrix(rnorm(new.samplesize, 0, error.sd), ncol=1)
  New.Ry <- matrix(mu, nrow = new.samplesize, ncol = 1) +
    New.Rxz %*% matrix(c(theta.x, rep(0, j.z)), ncol = 1) + Error
  New.OBS <- data.frame(New.Rxz[, col.index])
  colnames(New.OBS) <- c(paste('x', 1:j.x, sep=''), paste('z', 1:j.z, sep=''))[col.index]
  rm(Error)
  
  return(list(Rxz = Rxz, Ry = Ry, OBS = OBS, New.Rxz = New.Rxz, New.Ry = New.Ry, New.OBS = New.OBS))
}

Sim.Iter <- function(estimators = Estimators, samplesize = SampleSize, covariance = Covariance, error.sd = Error.SD, mu = Mu, theta.x =
                      Theta.X, j.x = J.X, j.z = J.Z, col.index = Col.Index, new.samplesize = New.SampleSize, include.intercept =
                      Include.Intercept, s.matrix = S.matrix, correct.model.form = Correct.Model.Form, potential.gamma =
                      Potential.Gamma, poppars = popPars, fix.x = Fix.X, fix.z = Fix.Z, potential.x = Potential.X, potential.z =
                      Potential.Z) {
  #####################################################################################
  #####################################################################################
  ########################
  ######################## INITIAL PART
  ########################
  #####################################################################################
  #####################################################################################
  
  #### I haven't set the random seed yet
  ## CHANGE THIS PART TO USE mvrnormArma
  nestimators <- length(estimators)
  Results <- data.frame(Bias = rep(NA, nestimators), MSE = rep(NA, nestimators), row.names = estimators)
  
  attach(poppars)
  randomData <- Data.Generation(samplesize = SampleSize, covariance = Covariance, error.sd = Error.SD, mu = Mu, theta.x = Theta.X, 
                               col.index = Col.Index, j.x = J.X, j.z = J.Z, new.samplesize = New.SampleSize) 
  attach(randomData)
  #####################################################################################
  #####################################################################################
  ########################
  ######################## MODEL AVERAGING PART
  ########################
  #####################################################################################
  #####################################################################################
  ## PREPARE ALL MODELS, OBTAIN A MATRIX WITH COLUMNS BEING ALL MODEL COEFFICIENT VECTORS
  OBS.y <- OBS[, 1]
  if (include.intercept == 1) {
    OBS.X <- cbind(1, OBS[, -1])
  } else {
    OBS.X <- OBS[, -1]
  }
  
  allModels <- EstAllModels(X = as.matrix(OBS.X), Xnew = as.matrix(cbind(1, New.OBS)), y = as.matrix(OBS.y), s = s.matrix)
  
  ## COMPUTE WEIGHTS FOR ALL METHODS
  AICw <- FMA(weighting.method = "AIC", allModels = allMods)
  BICw <- FMA(weighting.method = "BIC", allModels = allMods)
  JMAw <- FMA(weighting.method = "JMA", allModels = allMods)
  MMAw <- FMA(weighting.method = "MMA", allModels = allMods)
  
  ## COMPUTE PREDICTIONS
  Create.FMA(AICw, "AICmod", allModels)
  Create.FMA(BICw, "BICmod", allModels)
  Create.FMA(JMAw, "JMAmod", allModels)
  Create.FMA(MMAw, "MMAmod", allModels)
  
  ## COMPUTE BIAS AND MSE
  Results["SmoothAIC", "Bias"] <- mean(AICmod$preds - New.Ry)
  Results["SmoothAIC", "MSE"]  <- mean((AICmod$preds - New.Ry)^2)
  Results["SmoothBIC", "Bias"] <- mean(BICmod$preds - New.Ry)
  Results["SmoothBIC", "MSE"]  <- mean((BICmod$preds - New.Ry)^2)
  Results["JMA", "Bias"] <- mean(JMAmod$preds - New.Ry)
  Results["JMA", "MSE"]  <- mean((JMAmod$preds - New.Ry)^2)
  Results["MMA", "Bias"] <- mean(MMAmod$preds - New.Ry)
  Results["MMA", "MSE"]  <- mean((MMAmod$preds - New.Ry)^2)
  
  #### Ideal case where we know which regressors are important
  ######## Only meaningful when only finite number of regressors are used in the DGP. Not always valid!!
  # This is when we know the true model, we just don't know the coefficients
  if(samplesize >= (j.x + include.intercept)){ # If n >= p
    Results["Ideal", "Bias"] <-  mean(allModels$preds[, correct.model.form] - New.Ry)
    Results["Ideal", "MSE"] <-  mean((allModels$preds[, correct.model.form] - New.Ry)^2)
  }
  
  
  #### Model selection: The smallest AIC
  ModelChoose <- which.min(allModels$AIC)
  Results["Prac", "Bias"] <-  mean(allModels$preds[, ModelChoose] - New.Ry)
  Results["Prac", "MSE"] <-  mean((allModels$preds[, ModelChoose] - New.Ry)^2)
  
  #### Model selection: The largest possible model given the information at hand. Needed for adaptive penalization
  # POSSIBLE ERROR HERE, compare with old Pena.Weights[1:(length(fix.x)+length(fix.z))] <- 0 if fix.x=fix.z=NULL
  if (include.intercept == 1) {
    Adap.Pena.Weights <- 1/(allModels$coefficients[-1, M]^2)
    FullPreds <- as.matrix(cbind(1, New.OBS[, which(s.matrix[M, -1] == 1)])) %*% allModels$coefficients[, M]
  } else {
    Adap.Pena.Weights <- 1/(allModels$coefficients[, M]^2)
    FullPreds <- as.matrix(New.OBS[, which(s.matrix[M, ] == 1)]) %*% allModels$coefficients[, M]
  }
  
  Pena.Weights <- rep(1,length(fix.x)+length(potential.x)+length(fix.z)+length(potential.z))
  Pena.Weights[1:(length(fix.x)+length(fix.z))] <- 0
  Results["Full", "Bias"] <- mean(FullPreds - New.Ry)
  Results["Full", "MSE"] <- mean((FullPreds - New.Ry)^2)
  
  #####################################################################################
  #####################################################################################
  ########################
  ######################## PENALIZED PART
  ########################
  #####################################################################################
  #####################################################################################
  
  #### Lasso
  Results["Lasso", ] <- Wrapper.Glmnet(alphaval = 1,   penalty.weights = Pena.Weights,      methodname = "Lasso", include.intercept = include.intercept,
                                       ry = Ry, rxz = Rxz, col.index = col.index, new.OBS = New.OBS, new.Ry = New.Ry)
  Results["ALasso", ] <- Wrapper.Glmnet(alphaval = 1,   penalty.weights = Adap.Pena.Weights, methodname = "ALasso", include.intercept = include.intercept,
                                        ry = Ry, rxz = Rxz, col.index = col.index, new.OBS = New.OBS, new.Ry = New.Ry)
  Results["Ridge", ] <- Wrapper.Glmnet(alphaval = 0,   penalty.weights = Pena.Weights,      methodname = "Ridge", include.intercept = include.intercept,
                                       ry = Ry, rxz = Rxz, col.index = col.index, new.OBS = New.OBS, new.Ry = New.Ry)
  Results["ARidge", ] <- Wrapper.Glmnet(alphaval = 0,   penalty.weights = Adap.Pena.Weights, methodname = "ARidge", include.intercept = include.intercept,
                                        ry = Ry, rxz = Rxz, col.index = col.index, new.OBS = New.OBS, new.Ry = New.Ry)
  Results["ElasticNet", ] <- Wrapper.Glmnet(alphaval = 1/6, penalty.weights = Pena.Weights,      methodname = "ElasticNet", include.intercept = include.intercept,
                                            ry = Ry, rxz = Rxz, col.index = col.index, new.OBS = New.OBS, new.Ry = New.Ry)
  Results["AElasticNet", ] <- Wrapper.Glmnet(alphaval = 1/6, penalty.weights = Adap.Pena.Weights, methodname = "AElasticNet", include.intercept = include.intercept,
                                             ry = Ry, rxz = Rxz, col.index = col.index, new.OBS = New.OBS, new.Ry = New.Ry)
  
  #### SCAD
  SCAD.cv <- cv.ncvreg(X = OBS[, -1], y = Ry, family = c("gaussian"), penalty = "SCAD", gamma = 3.7, nlambda=100, eps=.001,
                       warn = FALSE, nfolds = 10, trace = FALSE, penalty.factor = Pena.Weights)
  SCAD.Pred <- predict(SCAD.cv, X = as.matrix(New.OBS), type='response', lambda = SCAD.cv$lambda.min)
  Results["SCAD", "Bias"] <-  mean(SCAD.Pred - New.Ry)
  Results["SCAD", "MSE"] <-  mean((SCAD.Pred - New.Ry)^2)
  rm(SCAD.cv)
  
  #### MCP
  for(i in 1:length(potential.gamma)){
    if(i == 1) {
      MCP.cv <- cv.ncvreg(X = OBS[, -1], y = Ry, family = c("gaussian"), penalty = "MCP", gamma = potential.gamma[i], nlambda=100,
                          eps = .001, warn = FALSE, nfolds = 10, trace = FALSE, penalty.factor = Pena.Weights)
      MCP.Min.CV <- MCP.cv$cve[which.min(MCP.cv$lambda)]
    } else {
      MCP.cv.New <- cv.ncvreg(X = OBS[,-1], y = Ry, family = c("gaussian"), penalty = "MCP", gamma = potential.gamma[i], nlambda = 100,
                              eps = .001, warn = FALSE, nfolds = 10, trace = FALSE, penalty.factor = Pena.Weights)
      if(MCP.Min.CV > MCP.cv.New$cve[which.min(MCP.cv.New$lambda)]){
        MCP.Min.CV <- MCP.cv.New$cve[which.min(MCP.cv.New$lambda)]
        MCP.cv <- MCP.cv.New
      }
    }
  }
  MCP.Pred <- predict(MCP.cv, X = as.matrix(New.OBS), type='response', lambda = MCP.cv$lambda.min)
  Results["MCP", "Bias"] <-  mean(MCP.Pred - New.Ry)
  Results["MCP", "MSE"]  <-  mean((MCP.Pred - New.Ry)^2)
  rm(MCP.cv)
  detach(poppars)
  detach(randomData)
  return(Results)
}

