\name{getFpTransforms}
\alias{getFpTransforms}
\title{Get the FP transforms matrix of a given covariate vector...}
\usage{getFpTransforms(vec, powers, center=TRUE)}
\description{Get the FP transforms matrix of a given covariate vector}
\details{Get the (centered) FP transforms matrix of a given covariate vector, when the corresponding 
power vector (with at least one element) is provided.}
\value{the FP transforms matrix with proper colnames.}
\keyword{utilities}
\keyword{internal}
\arguments{\item{vec}{positive (== already shifted and scaled) column vector with proper colname}
\item{powers}{power vector with at least one element}
\item{center}{center the columns of the FP transforms matrix around zero? (default)}}
