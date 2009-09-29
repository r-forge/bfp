\name{scrHpd}
\alias{scrHpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~functions to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the functions
     scrHpd do. ~~
}
\usage{
scrHpd(samples, mode = apply(samples, 2, median), level = 0.95, grid = ncol(samples))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{ ~~Describe \code{samples} here~~ }
  \item{mode}{ ~~Describe \code{mode} here~~ }
  \item{level}{ ~~Describe \code{level} here~~ }
  \item{grid}{ ~~Describe \code{grid} here~~ }
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
## scrHpd
function(samples,          # sample matrix
                   mode = apply (samples, 2, median), # mode vector of length = ncol (samples)
                   level = 0.95,     # credibility level
                   grid = ncol (samples) # number of elements in each sample vector
                   )
{
    n <- ncol(samples)
    m <- nrow(samples)

    if (n != length (mode))
        stop ("mode vector must have same length as samples matrix!")

    distance <- abs (sweep (samples, 2, mode)) # absolute distance from mode vector

    k <- floor (level * m)
    ## Calculates a simultaneous (k/m)*100\% credible band
    ## using the ranks approach

    rankdistance <- apply (distance, 2, rank) # colwise (= elementwise) ranks of distances

    tstari <- apply (rankdistance, 1, max) # maximum ranks in each sample
    ordtstari <- sort.int (tstari, method = "quick") # and sort the maximum ranks

    tstar <- ordtstari[k]               # the required rank

    selectMat <- rankdistance <= tstar

    ret <- matrix (nrow = 2, ncol = n)
    for (i in seq_len (n)){
        ret[, i] <- range (samples[selectMat[,i], i])
    }

    return (ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{internal}
