\name{Wrapper.Glmnet}
\alias{Wrapper.Glmnet}
\title{Wrapper glmnet}
\description{
  Wrapper for \code{\link{glmnet}}.
}
\usage{
Wrapper.Glmnet(alphaval, penalty.weights, methodname, 
include.intercept = Include.Intercept, ry = Ry, rxz = Rxz, col.index = Col.Index, 
new.OBS = New.OBS, new.Ry = New.Ry)
}
\arguments{
  \item{alphaval}{Value for alpha in glmnet.}
  \item{penalty.weights}{Value for penalty weights in glmnet.}
  \item{methodname}{Name of the method.}
  \item{include.intercept}{If intercept is included.}
  \item{ry}{The response variable.}
  \item{rxz}{The regressor matrix.}
  \item{col.index}{Column index.}
  \item{new.OBS}{Matrix of new observations.}
  \item{new.Ry}{New responses.}
}
\keyword{Wrapper.Glmnet}