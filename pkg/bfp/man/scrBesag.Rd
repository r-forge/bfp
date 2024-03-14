\name{scrBesag}
\alias{scrBesag}
\title{Simultaneous credible band computation (Besag, Green et al algorithm)}
\usage{scrBesag(samples, level=0.95)}
\description{Simultaneous credible band computation}
\details{Calculates a series of simultaneous credible bounds for one parameter type, following section 6.3 in the seminal paper "Bayesian computation and stochastic systems". The corresponding algorithm was invented by Peter Green.}
\value{matrix with \sQuote{lower} and \sQuote{upper} rows}
\references{J. Besag, P. Green, D. Higdon, K. Mengersen (1995):
Bayesian computation and stochastic systems,
\emph{Statistical Science} \bold{10}/1, 3--41}
\author{Thomas Kneib}
\keyword{htest}
\arguments{\item{samples}{m by n matrix where m is the number of parameters,
n is the number of samples and  hence each (multivariate) sample is a column
in the matrix \code{samples}}
\item{level}{the credible level (default: 0.95)}}
