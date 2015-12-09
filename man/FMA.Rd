\name{FMA}
\alias{FMA}
\title{Frequentist Model Averaging}
\description{
  Frequentist model averaging using smoothed AIC, smoothed BIC, Mallow's model averaging or Jackknife model averaging.
}
\usage{
FMA(weighting.method, allMods, X, y, submodels, include.intercept = 1)
}
\arguments{
  \item{weighting.method}{One of \code{"AIC"} (for smoothed AIC), \code{"BIC"} (for smoothed BIC), \code{"MMA"} (for Mallow's model averaging) and \code{"JMA"} (for Jackknife model averaging).}
  \item{allMods}{The list returned by \code{EstAllModels}.}
  \item{X}{Matrix of independent variables.}
  \item{y}{Vector of the dependent variable.}
  \item{submodels}{Either one of \code{"all"}, \code{"singleton"} or \code{"nested"} (which is then passed to \code{Generate.S}), or a custom selection matrix (of appropriate size).}
  \item{include.intercept}{Flag for inclusion of intercept. Default is 1 for inclusion.}
}
\details{Either one must supply the function with the \code{allMods} argument, or the raw data and what models to average over. If one gives the function \code{allMods}, then only the weighting.method needs to be specified. If \code{allMods} is not supplied, then one must give the function the raw data and the selection matrix so that it can call \code{EstAllModels}. The first alternative is to be preferred if one will use multiple averaging techniques for the same data. See the examples.}
\value{Returns a vector with the weights.}
\seealso{
  \code{\link{EstAllModels}},
    \code{\link{Generate.S}}.
}
\examples{
## First alternative: supply allMods
\dontrun{myAllMods <- EstAllModels(X = Xdata, Xnew = Xnewdata, y = ydata, s = smatrix)}
\dontrun{wAIC <- FMA("AIC", allMods = myAllMods)}
\dontrun{wJMA <- FMA("JMA", allMods = myAllMods)}

## Second alternative: let the function call EstAllModels
\dontrun{wAIC <- FMA("AIC", X = Xdata, Xnew = Xnewdata, y = ydata, submodels = "nested")}
\dontrun{wJMA <- FMA("JMA", X = Xdata, Xnew = Xnewdata, y = ydata, submodels = "nested")}
}
\keyword{FMA}