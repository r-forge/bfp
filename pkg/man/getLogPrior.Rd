\name{getLogPrior}
\alias{getLogPrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~functions to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the functions
     getLogPrior do. ~~
}
\usage{
getLogPrior(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
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
## getLogPrior
function (x       # a valid BayesMfp-Object of length 1 (otherwise only first element recognized)
                         )
{
    x <- x[1]
    powers <- x[[1]]$powers
    maxs <- attr (x, "shiftScaleMax")[, 3]
    cards <- attr (x, "shiftScaleMax")[, 4]

    logVals <- numeric (0)
    for (i in seq_along (powers)){
        deg <- length (powers[[i]])
        logVals[i] <- - lchoose(cards[i] - 1 + deg, deg) - log1p(maxs[i])
    }

    return (sum (logVals) - max (attr (x, "indices")$uc) * log (2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{internal}

