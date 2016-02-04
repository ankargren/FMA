\name{PIE}
\alias{PIE}
\title{Perturbation index}
\description{
  Perturbation index (Yuan and Yang (2005)) CURRENTLY DOES NOT WORK
}
\usage{
PIE(y, X, rhovec = seq(0.2, 0.8, by = 0.2), nreps = 100, avertype = "all", 
    estimators = c("AIC", "BIC", "MMA", "JMA"))
}
\arguments{
  \item{y}{Response vector.}
  \item{X}{Design matrix.}
  \item{rhovec}{Vector of rho values.}
  \item{nreps}{Number of replications.}
  \item{avertype}{Type of averaging.}
  \item{estimators}{The (FMA) estimators to include.}
}
\value{Returns matrix with the index}
\keyword{PIE}