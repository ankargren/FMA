Wrapper.Glmnet <- function(alphaval, Xdata, ydata, Xnewdata, mu) {
  cv <- cv.glmnet(x = Xdata, y = ydata, nfolds = 10, alpha = alphaval, nlambda = 100, type.measure = "mse",
                  standardize = TRUE, intercept = 1)
  Pred <- predict(cv, newx = Xnewdata, s = "lambda.min")
  return(mean((mu - Pred)^2))
}