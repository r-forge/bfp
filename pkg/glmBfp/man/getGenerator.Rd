\name{getGenerator}
\alias{getGenerator}
\title{Internal helper function which gets the generator (and normalizing...}
\usage{getGenerator(method, unnormDensFun, xVals, unnormDensVals)}
\description{Internal helper function which gets the generator (and normalizing
constant)}
\value{a list with the elements \dQuote{generator} and \dQuote{normConst},
containing the generator function (with argument \code{n} for the number
of samples) and the normalizing constant of the density, respectively.}
\keyword{internal}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{method}{the methods string}
\item{unnormDensFun}{the unnormalized density function}
\item{xVals}{the x values}
\item{unnormDensVals}{the unnormalized density values}}
