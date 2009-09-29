\name{inclusionProbs}
\alias{inclusionProbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inclusionProbs(BayesMfpObject, postProbs = posteriors(BayesMfpObject, ind = 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BayesMfpObject}{ ~~Describe \code{BayesMfpObject} here~~ }
  \item{postProbs}{ ~~Describe \code{postProbs} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (            # compute posterior inclusion probabilites based on BayesMfpObject
                            BayesMfpObject,
                            postProbs = posteriors (BayesMfpObject, ind = 1)
                            )
{
    postProbs <- postProbs / sum (postProbs)
    inds <- attr (BayesMfpObject, "indices")
    termNames <- attr (BayesMfpObject, "termNames")

    nams <- unlist (termNames[c ("bfp", "uc")])
    ret <- numeric (length (nams))
    names (ret) <- nams

    i <- 0

    for (j in seq_along (inds$bfp)){
        i <- i + 1
        present <- sapply (BayesMfpObject, function (one) as.logical (length (one$powers[[j]])))
        ret[i] <- sum (present * postProbs)
    }

    for (j in seq_along (inds$ucList)){
        i <- i + 1
        present <- sapply (BayesMfpObject, function (one) any (j == one$ucTerms))
        ret[i] <- sum (present * postProbs)
    }

    return (ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
