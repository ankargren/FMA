\name{Hansen}
\alias{Hansen}
\title{Generate data as in Hansen (2007)}
\description{
  This function generates data as in Hansen (2007).
}
\usage{
Hansen(alpha, n, r2, J)
}
\arguments{
  \item{alpha}{Scalar determining the rate of decay.}
  \item{n}{The sample size.}
  \item{r2}{The population coefficient of determination.}
  \item{J}{The number of variables to use to approximate the infinite sum.}
}
\value{
Returns a list with the following:
\item{X}{An \eqn{n * J} matrix containing the independent variables.}
\item{y}{An \eqn{n * 1} matrix containing the dependent variable.}
\item{mu}{An \eqn{n * 1} matrix containing the true conditional mean.}
}
\examples{
## Generate data as in Hansen:
\notrun{data <- Hansen(alpha = 0.5, n = 50, r2 = 0.5, J = 1000)}
}
\references{
Hansen, Bruce (2007) Least Squares Model Averaging, Econometrica, 75(4):1175-1189.
}
\keyword{Hansen}