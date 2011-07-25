\name{cppOptimize}
\alias{cppOptimize}
\title{Interface to the internal C++ optimization routine "optimize"...}
\usage{cppOptimize(f_, min.x, max.x, prec=.Machine$double.eps^0.25)}
\description{Interface to the internal C++ optimization routine "optimize"}
\value{A list with the following elements:
\describe{
\item{par}{the minimum abscissa found by the algorithm}
\item{inv.hessian}{the inverse Hessian at \code{par}}
\item{evaluations}{list of the function evaluation pairs: \code{args} and
\code{vals}} 
}}
\keyword{utilities}
\keyword{internal}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{f_}{the function}
\item{min.x}{minimum bound on x}
\item{max.x}{maximum bound on x}
\item{prec}{precision (same default as for the original optimize function)}}
