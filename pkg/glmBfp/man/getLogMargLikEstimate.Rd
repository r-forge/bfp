\name{getLogMargLikEstimate}
\alias{getLogMargLikEstimate}
\title{Compute the Chib-Jeliazkov log marginal likelihood estimate from the MCMC output...}
\usage{getLogMargLikEstimate(numeratorTerms, denominatorTerms,
    highDensityPointLogUnPosterior, endLag=40)}
\description{Compute the Chib-Jeliazkov log marginal likelihood estimate from the MCMC output}
\value{list with resulting \dQuote{estimate} and \dQuote{standardError}.}
\keyword{internal}
\arguments{\item{numeratorTerms}{terms for the sum in the numerator of the posterior density ordinate
estimate at the high density point}
\item{denominatorTerms}{terms for the sum in the denominator of the posterior density ordinate 
estimate at the high density point}
\item{highDensityPointLogUnPosterior}{log unnormalized posterior (i.e. likelihood times prior)
of the high density point}
\item{endLag}{up to which lag should the covariance incorporate the uncertainty? (default: 40)}}
