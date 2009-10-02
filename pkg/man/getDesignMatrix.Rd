\name{getDesignMatrix}
\alias{getDesignMatrix}

\title{Extract the design matrix of a multiple FP model.}
\description{
  Extract the (centered) design matrix of the first element (== model) of a
  \code{\link{BayesMfp}} object. This is an internal function not
  intended to be used publicly.
}
\usage{
getDesignMatrix(x, fixedColumns=TRUE)
}
\arguments{
  \item{x}{a valid BayesMfp-Object of length 1 (otherwise only first
    element is recognized)}
  \item{fixedColumns}{return the fixed columns (intercept) inside the
    matrix (default) or not?}  
}
\value{
  The design matrix, where the non-fixed part is columnwise centered
  (that is, the colmeans are zero).
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{getFpTransforms}}}
\examples{
## construct a BayesMfp object
set.seed(19)

x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)

y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 200, method="exhaustive")

## get the design matrix of the fifth best model
bfp:::getDesignMatrix(test[5])
}
\keyword{utilities}
\keyword{internal}
