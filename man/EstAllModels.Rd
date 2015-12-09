\name{EstAllModels}
\alias{EstAllModels}
\title{Estimate all submodels}
\description{
  This function estimates all submodels.
}
\usage{
EstAllModels(X, Xnew, y, s)
}
\arguments{
  \item{X}{A design matrix.}
  \item{Xnew}{Design matrix used for predictions.}
  \item{y}{Vector of the response.}
  \item{s}{Matrix where the rows correspond to a submodel.}
}
\value{
Returns a list with the following:
\item{SEs}{Variances of the submodels.}
\item{coefficients}{Coefficients from each submodel.}
\item{K}{The number of included variables in each submodel.}
\item{etilde}{Matrix of residuals (for MMA). Each column corresponds to residuals from a submodel.}
\item{eJMA}{The matrix of error terms for JMA. Each column corresponds to jackknife residuals in the submodels.}
\item{AIC}{The AIC values for the models.}
\item{BIC}{The BIC values for the models.}
\item{preds}{The predictions from all submodels.}
}
\seealso{
    \code{\link{Generate.S}}.
}
\examples{
## Generate some data
Xdata <- matrix(rnorm(100), 25, 4)
ydata <- rowSums(Xdata) + rnorm(25)
Xnewdata <- matrix(rnorm(100), 25, 4)

## Get the submodel matrix
smatrix <- Generate.S(Fix.X = 1, Potential.X = 3, avertype = "nested")

## Estimate all submodels
allMods <- EstAllModels(X = Xdata, Xnew = Xnewdata, y = ydata, s = smatrix)
}
\keyword{EstAllModels}