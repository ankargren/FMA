\name{Hansen}
\alias{Hansen}
\title{Generate data as in Hansen (2007)}
\description{
  This function generates data as in Hansen (2007).
}
\usage{
Hansen(alpha, n, r2, J, JJ = 1000, rho = 0, equi = FALSE, include.intercept = TRUE)
}
\arguments{
  \item{alpha}{Scalar determining the rate of decay.}
  \item{n}{The sample size.}
  \item{r2}{The population coefficient of determination.}
  \item{J}{The number of variables to observe.}
  \item{JJ}{The number of terms to use to approximate infinite sums.}
  \item{rho}{The correlation.}
  \item{equi}{Logical, if \code{TRUE} then the covariance matrix is an equicorrelation matrix with correlation coefficient \code{rho}. If \code{FALSE}, then the correlation is power decaying.}
  \item{include.intercept}{Logical, if \code{TRUE} then an intercept is included and there are \code{J - 1} random variables.}
}
\value{
Returns a list with the following:
\item{X}{An \eqn{n * J} matrix containing the independent variables.}
\item{y}{An \eqn{n * 1} matrix containing the dependent variable.}
\item{mu}{An \eqn{n * 1} matrix containing the true conditional mean.}
}
\examples{
## Generate data as in Hansen:
\dontrun{data <- Hansen(alpha = 0.5, n = 50, r2 = 0.5, J = ceiling(3*50^(1/3)), JJ = 1000, rho = 0, equi = FALSE)}
}
\references{
Hansen, Bruce (2007) Least Squares Model Averaging, Econometrica, 75(4):1175-1189.
}
\keyword{Hansen}