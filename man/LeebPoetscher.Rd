\name{LeebPoetscher}
\alias{LeebPoetscher}
\title{Generate data as in Leeb and Pötscher (2006)}
\description{
  Data generation corresponding to Setups I-III by Leeb and Pötscher (2006).
}
\usage{
LeebPoetscher(samplesize, gam, eta, setup)
}
\arguments{
  \item{samplesize}{The sample size.}
  \item{gam}{The drifting parameter.}
  \item{setup}{Either 1, 2 och 3 for the setup to use.}
}
\value{
Returns a list with the following:
\item{X}{An \eqn{n * 8} matrix containing the independent variables}
\item{y}{An \eqn{n * 1} matrix containing the dependent variable.}
\item{mu}{An \eqn{n * 1} matrix containing the true conditional mean.}
}
\keyword{LeebPoetscher}