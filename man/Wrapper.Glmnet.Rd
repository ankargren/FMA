\name{Wrapper.Glmnet}
\alias{Wrapper.Glmnet}
\title{Wrapper function for using glmnet}
\description{
  This function is a wrapper for choosing the bandwidth in glmnet and then predicting using this choice.
}
\usage{
Wrapper.Glmnet(alphaval, Xdata, ydata, Xnewdata, mu) 
}
\arguments{
  \item{alphaval}{The choice of \code{alpha} in \code{\link{glmnet}}. 0 corresponds to the Lasso, 1 to ridge.}
  \item{Xdata}{The design matrix. Inclusion of a constant vector is not necessary.}
  \item{ydata}{The response vector.}
  \item{Xnewdata}{The data to use for prediction.}
  \item{mu}{The true mean vector to use for computing MSE of predictions.}
}
\value{
Returns a scalar value with the MSE of the prediction of the mean vector. 
}
\keyword{Wrapper.Glmnet}
\seealso{
\code{\link{glmnet}}, \code{\link{cv.glmnet}}.
}