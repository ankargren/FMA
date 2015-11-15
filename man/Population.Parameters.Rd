\name{Population.Parameters}
\alias{Population.Parameters}
\title{Create the population parameters}
\description{
  Creates the population parameters as in Hansen (2007) given the inputs.
}
\usage{
Population.Parameters(J.X, J.Z, Cov.Const, Var.Const, Include.Intercept)
}
\arguments{
  \item{J.X}{Number of Xs.}
  \item{J.Z}{Number of Zs.}
  \item{Cov.Const}{A constant which is used in the correlation matrix. At most, the correlation is this number.}
  \item{Var.Const}{The variance of the independent variables.}
}
\keyword{Population.Parameters}