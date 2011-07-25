#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[GlmBayesMfp-methods.R] by DSB Die 25/05/2010 16:08 (CEST)>
##
## Description:
## Additional convenience methods for GlmBayesMfp class objects.
##
## History:
## 12/02/2010   modify from bfp package
## 15/03/2010   adapt print.GlmBayesMfp to new S4 gPrior element
#####################################################################################

##' @include posteriors.R
{}

##' Extract method for GlmBayesMfp objects
##'
##' Extract a subset of models from a \code{\link{GlmBayesMfp}} object.
##'
##' @S3method "[" GlmBayesMfp
##' @method [ GlmBayesMfp
##' @param x valid \code{\link{GlmBayesMfp}} object
##' @param \dots transports the indexes of the models
##' @return The subsetted object.
##'
##' @name Extract.GlmBayesMfp
##' @aliases Extract.GlmBayesMfp [.GlmBayesMfp
##' @keywords methods
##' @seealso \code{\link{glmBayesMfp}} 
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
'[.GlmBayesMfp' <- function (x, ...)       
{
    y <- NextMethod("[")
    mostattributes (y) <- attributes (x)
    names (y) <- names (x)[...]
    class(y) <- oldClass(x)
    y
}



##' Print a GlmBayesMfp object.
##'
##' @S3method print GlmBayesMfp
##' @method print GlmBayesMfp
##' @param x valid \code{\link{GlmBayesMfp}} object
##' @param \dots unused
##' @return Only used for its side effect
##' 
##' @keywords methods
##' @seealso \code{\link{glmBayesMfp}} 
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
print.GlmBayesMfp <- function (x, ...)
{
    cat ("------------------------------ GlmBayesMfp-Output ------------------------------\n")
    cat (length (x), "multivariable fractional polynomial model(s) of total (visited/cached)",
         attr (x, "numVisited"), "for following covariates:\n\n")
    cat ("fixed:                ", paste(attr (x, "termNames")$fixed, collapse = ", "), "\n")
    cat ("uncertain fixed form: ", paste(attr (x, "termNames")$uc, collapse = ", "), "\n")
    cat ("fractional polynomial:", paste(attr (x, "termNames")$bfp, collapse = ", "), "\n")
    pr <- attr (x, "priorSpecs")
    cat ("\nPrior for g was,", class(pr$gPrior),
         "\nand a ", pr$modelPrior, " prior on the model space was used.",
         "\n")
}


##' Convert a GlmBayesMfp object into a data frame
##'
##' @S3method as.data.frame GlmBayesMfp
##' @method as.data.frame GlmBayesMfp
##' @param x valid \code{\link{GlmBayesMfp}} object
##' @param row.names optional rownames (default is to keep the names of the
##' \code{\link{GlmBayesMfp}} list) 
##' @param \dots unused
##' @param freq should empirical frequencies of the models in the sampling
##' path be given? (default)
##' @return The data frame with the following columns:
##' \describe{
##'     \item{posterior}{the posterior model probabilities}
##'     \item{logMargLik}{the log marginal likelihood of the models}
##'     \item{logPrior}{the log prior probabilities of the models}
##' }
##' Additionally, for each uncertain fixed form covariates a column with the inclusion
##' status, and for each fractional polynomial a column with the powers are returned.
##' @seealso \code{\link{glmBayesMfp}} 
##' @keywords methods
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
as.data.frame.GlmBayesMfp <- function (x, 
                                       row.names = NULL,
                                       ...,        
                                       freq = TRUE 
                                       )
{
    ## posterior probabilites:
    ret <- data.frame (posterior=
                       posteriors(x, type="normalized"))
    if (! is.null(attr (x, "chainlength")) &&
        freq)
        ret$frequency <- posteriors (x, type="sampling")

    ## row names
    row.names (ret) <-
        if(! is.null(row.names))
            row.names
        else
            names(x)    

    ## log marginal likelihood and log prior
    ret$logMargLik <- logMargLiks(x)
    ret$logPrior <- logPriors(x)

    ## uncertain fixed form covariates:
    for (i in seq_along (ucNames <- attr (x, "termNames")$uc))
    {
        ret[[ucNames[i]]] <- sapply (x,
                                     function(one)
                                     i %in% one$configuration$ucTerms)
    }
    ## fractional polynomial:
    for (fpName in attr (x, "termNames")$bfp)
    {
        ret[[fpName]] <- sapply (x,
                                 function(one)
                                 paste(one$configuration$powers[[fpName]],
                                       collapse = ", "))
    }

    ret
}
