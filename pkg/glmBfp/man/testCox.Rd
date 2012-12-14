\name{testCox}
\alias{testCox}
\title{Test the Cox model computation for the TBF approach}
\usage{
  testCox(survTimes, censInd, X, useCppCode = FALSE)
}
\arguments{
  \item{survTimes}{the survival times}

  \item{censInd}{the logical censoring indicator
  (\code{TRUE} = observed, \code{FALSE} = censored survival
  times)}

  \item{X}{the design matrix, *without* the intercept 1's!,
  with at least one column}

  \item{useCppCode}{use the C++ code? (default) otherwise
  the R-function \code{\link{coxph}} is used}
}
\value{
  a list with the coefficient estimates (\code{betas}), the
  covariance matrix estimate (\code{cov}) and the residual
  deviance (\code{deviance}).
}
\description{
  Test the Cox model computation for the TBF approach
}
\author{
  Daniel Sabanes Bove
}
\keyword{internal}
\keyword{utilities}

