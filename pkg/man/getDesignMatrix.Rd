\name{getDesignMatrix}
\alias{getDesignMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~functions to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the functions
     getDesignMatrix do. ~~
}
\usage{
getDesignMatrix(x, fixedColumns=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a valid BayesMfp-Object of length 1 (otherwise only first
    element is recognized)}
  \item{fixedColumns}{return the fixed columns inside the matrix or
    not?} 
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  The design matrix.
}
\references{ ~put references to the literature/web site here ~ }
\author{Daniel Saban\'es Bov\'e}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The functions are currently defined as
## getDesignMatrix
function (x          # a valid BayesMfp-Object of length 1 (otherwise only first element recognized)
                             )
{
    full <- attr (x, "x")
    inds <- attr (x, "indices")
    powers <- x[1][[1]]$powers
    ucSet <- x[1][[1]]$ucTerms


    nFix <- length (inds$fixed)

    ucColInds <- inds$uc \%in\% ucSet
    nUc <- sum (ucColInds)

    nFp <- length (unlist (powers))

    ret <- matrix (nrow = nrow (full), ncol = nFix + nUc + nFp)
    retColnames <- character (ncol (ret))
    col <- 0                            # invariant: already col columns written

    ## fixed columns
    new <- full[, inds$fixed, drop = FALSE]
    newInds <- col + seq_along (inds$fixed)

    ret[, newInds] <- new
    retColnames[newInds] <- colnames (new)

    col <- col + nFix

    ## fp part
    for (i in seq_along (inds$bfp)){
        pi <- powers[[i]]
        if (len <- length (pi)) {       # if there is at least one power
            new <- getFpTransforms (full[, inds$bfp[i], drop = FALSE], pi)
            newInds <- col + seq_along (pi)

            ret[, newInds] <- new
            retColnames[newInds] <- colnames (new)

            col <- col + len
        }
    }

    ## uc part
    if (length (ucSet)){
        new <- full[, ucColInds, drop = FALSE]
        newInds <- col + seq_len (nUc)

        ret[, newInds] <- new
        retColnames[newInds] <- colnames (new)

        col <- col + nUc
    }

    rownames (ret) <- rownames (full)
    colnames (ret) <- retColnames

    return (ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{internal}
