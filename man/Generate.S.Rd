\name{Generate.S}
\alias{Generate.S}
\title{Generate the S matrix}
\description{
  Generation of the S matrix. 
}
\usage{
Generate.S(Fix.X, Potential.X, avertype = "all")
}
\arguments{
  \item{Fix.X}{Number of fixed Xs, i.e. variables which should always be included.}
  \item{Potential.X}{Number of potential variables, i.e. variables which may or may not be included.}
  \item{avertype}{Type of averaging to be done. Can be either \code{singleton} (models with only one regressor at a time), \code{nested} (all nested models) or \code{all} (all combinations of submodels).}
}
\value{Returns a selection matrix.}
\examples{
## all nested models with 5 variables, first variable always included
smatrix <- Generate.S(Fix.X = 1, Potential.X = 5, avertype = "nested")

## use this together with EstAllModels to estimate all the submodels efficiently
\dontrun{allMods <- EstAllModels(X = Xdata, Xnew = Xnewdata, y = ydata, s = smatrix)}
}
\keyword{Generate.S}