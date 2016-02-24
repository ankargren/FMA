Wrapper.Glmnet <- function(alphaval, Xdata, ydata, Xnewdata, mu, prediction = TRUE, include.intercept = 1) {
  cv <- cv.glmnet(x = Xdata, y = ydata, nfolds = 10, alpha = alphaval, nlambda = 100, type.measure = "mse",
                  standardize = TRUE, intercept = include.intercept)
  if (prediction == TRUE) {
    Pred <- predict(cv, newx = Xnewdata, s = "lambda.min")
    ret <- mean((mu - Pred)^2)
  } else {
    ret <- as.matrix(coef(cv, s = "lambda.min"))[-1]
  }
  return(ret)
}