\name{myrmvt}
\alias{myrmvt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~functions to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the functions
     myrmvt do. ~~
}
\usage{
myrmvt(n, sigma = diag(2), mu = rep(0, 2), df = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{sigma}{ ~~Describe \code{sigma} here~~ }
  \item{mu}{ ~~Describe \code{mu} here~~ }
  \item{df}{ ~~Describe \code{df} here~~ }
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
## myrmvt
function (n, sigma = diag (2), mu = rep (0, 2), df = 1)
{
    ret <- rmvt (n = n, sigma = sigma, df = df)
    ret <- sweep (ret, 2, mu, "+")
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{internal}% __ONLY ONE__ keyword per line
