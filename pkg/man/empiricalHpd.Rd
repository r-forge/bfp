\name{empiricalHpd}
\alias{empiricalHpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~functions to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the functions
     empiricalHpd do. ~~
}
\usage{
empiricalHpd(theta, level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{ ~~Describe \code{theta} here~~ }
  \item{level}{ ~~Describe \code{level} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the values returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{Daniel Saban\'es Bov\'e}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The functions are currently defined as
## empiricalHpd
function (theta,        # sample vector of parameter
                          level         # credibility level
                          )
{
    M <- length (theta)
    thetaorder <- sort.int (theta, method = "quick")

    alpha <- 1 - level
    maxSize <- ceiling (M * alpha)
    ind <- seq_len (maxSize)

    lower <- thetaorder[ind]
    upper <- thetaorder[M - maxSize + ind]
    size <- upper - lower
    sizeMin <- which.min(size)

    HPDLower <- thetaorder[sizeMin]
    HPDUpper <- thetaorder[M - maxSize + sizeMin]
    return (c (HPDLower, HPDUpper))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{internal}
