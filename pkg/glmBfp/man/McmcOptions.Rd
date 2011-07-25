\name{McmcOptions}
\alias{McmcOptions}
\title{Constructor for class McmcOptions...}
\usage{McmcOptions(iterations=as.integer(burnin + (step * samples)),
    burnin=10000, step=10, samples=10000)}
\description{Constructor for class McmcOptions}
\details{Note that the argument \code{samples} is included for convenience only -
you can specify it instead of \code{iterations}.}
\value{the freshly built object of class \code{\linkS4class{McmcOptions}}}
\keyword{classes}
\arguments{\item{iterations}{number of MCMC iterations (default: \code{110,000})}
\item{burnin}{number of burn-in iterations which are not saved (default:
\code{10,000})}
\item{step}{only every step-th iteration is saved after the burn-in
(default: \code{10})}
\item{samples}{number of resulting samples (by default \code{10,000} will result)}}
