\name{dinvGamma}
\alias{dinvGamma}
\alias{pinvGamma}
\alias{qinvGamma}
\alias{rinvGamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~functions to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the functions
     dinvGamma, pinvGamma, qinvGamma, rinvGamma do. ~~
}
\usage{
dinvGamma(x, a, b, log = FALSE, normalize = TRUE)
pinvGamma(q, a, b, lower.tail = TRUE, log.p = FALSE)
qinvGamma(p, a, b, lower.tail = TRUE, log.p = FALSE)
rinvGamma(n, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{a}{ ~~Describe \code{a} here~~ }
  \item{b}{ ~~Describe \code{b} here~~ }
  \item{log}{ ~~Describe \code{log} here~~ }
  \item{normalize}{ ~~Describe \code{normalize} here~~ }
  \item{q}{ ~~Describe \code{q} here~~ }
  \item{lower.tail}{ ~~Describe \code{lower.tail} here~~ }
  \item{log.p}{ ~~Describe \code{log.p} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{n}{ ~~Describe \code{n} here~~ }
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
## dinvGamma
function (x, a, b, log = FALSE, normalize = TRUE)
{
    ret <- - (a + 1) * log (x) - b / x
    if (normalize)
        ret <- ret + a * log (b) - lgamma (a)
    if (log)
        return (ret)
    else
        return (exp (ret))
  }
## pinvGamma
function (q, a, b, lower.tail = TRUE, log.p = FALSE)
{
    pgamma (q = 1 / q, shape = a, rate = b, lower.tail = !lower.tail, log.p = log.p)
  }
## qinvGamma
function (p, a, b, lower.tail = TRUE, log.p = FALSE)
{
    1 / qgamma (p = p, shape = a, rate = b, lower.tail = !lower.tail, log.p = log.p)
  }
## rinvGamma
function (n, a, b)
{
    1 / rgamma (n, shape = a, rate = b)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{internal}
