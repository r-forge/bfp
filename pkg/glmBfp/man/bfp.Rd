\name{bfp}
\alias{bfp}
\title{Mark a covariate for transformation with fractional polynomials...}
\usage{bfp(x, max=2, scale=TRUE, rangeVals)}
\description{Mark a covariate for transformation with fractional polynomials}
\details{Using this function, you can mark a covariate for transformation with fractional polynomials
in the formula for \code{\link{glmBayesMfp}}.}
\value{name of the provided covariate, with the other function parameters as  attached attributes}
\keyword{utilities}
\arguments{\item{x}{the covariate name}
\item{max}{maximum degree for this FP}
\item{scale}{use pre-transformation shifting and scaling to avoid numerical problems?}
\item{rangeVals}{extra numbers if the shifting and scaling should consider values in this range}}
