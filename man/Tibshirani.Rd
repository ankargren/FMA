\name{Tibshirani}
\alias{Tibshirani}
\title{Generate data as in Tibshirani (1996)}
\description{
  Data generation corresponding to the four simulation setups used by Tibshirani (1996).
}
\usage{
Tibshirani(samplesize, ex, equi = FALSE)
}
\arguments{
  \item{samplesize}{The sample size.}
  \item{ex}{The example to generate data according to. Can be 1, 2, 3 or 4.}
  \item{equi}{True or false to indicate if equicorrelations should be used.}
}
\value{
Returns a list with the following:
\item{X}{An \eqn{n * p} matrix containing the independent variables (\code{p} is 8 for \code{ex \%in\% 1:3}, 40 for \code{ex == 4}.}
\item{y}{An \eqn{n * 1} matrix containing the dependent variable.}
\item{mu}{An \eqn{n * 1} matrix containing the true conditional mean.}
}
\keyword{Tibshirani}