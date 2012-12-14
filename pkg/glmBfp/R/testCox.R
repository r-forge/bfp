#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian fractional polynomials
##        
## Time-stamp: <[testCox.R] by DSB Mon 03/12/2012 18:21 (CET)>
##
## Description:
## R code and interface to C++ code that could be used for the Cox model
## computations, for regression testing purposes.
##
## History:
## 29/11/2012   file creation
## 03/12/2012   make it internal
#####################################################################################

##' @include helpers.R
{}


##' Test the Cox model computation for the TBF approach
##'
##' @param survTimes the survival times 
##' @param censInd the logical censoring indicator (\code{TRUE} = observed, 
##' \code{FALSE} = censored survival times) 
##' @param X the design matrix, *without* the intercept 1's!, with at least 
##' one column
##' @param useCppCode use the C++ code? (default) otherwise the R-function
##' \code{\link{coxph}} is used
##' @return a list with the coefficient estimates (\code{betas}), the
##' covariance matrix estimate (\code{cov}) and the residual deviance
##' (\code{deviance}).
##'
##' @keywords utilities internal
##' @importFrom survival coxph
##' @author Daniel Sabanes Bove
testCox <- function(survTimes,
                    censInd,
                    X,
                    useCppCode=FALSE)
{
    ## checks
    stopifnot(is.numeric(survTimes),
              is.logical(censInd),
              is.matrix(X),
              identical(length(survTimes), length(censInd)),
              identical(length(survTimes), nrow(X)),
              is.bool(useCppCode))

    ## summaries
    nObs <- nrow(X)
    nCovs <- ncol(X)

    ## check that we have at least one covariate
    ## (null model is not treated here)
    stopifnot(nCovs > 0)

    ## method for handling ties:
    method <- "efron"
    
    ## now do either simple R or go C++
    if(useCppCode)
    {
        ## sort according to survival times
        sorted <- order(survTimes)
        survTimes <- survTimes[sorted]
        censInd <- censInd[sorted]
        X <- X[sorted, ]

        ## go C++ for fitting the Cox model
        storage.mode(X) <- "double"
        tmp <- .External(cpp_coxfit,
                         as.double(survTimes),
                         as.integer(censInd),
                         X,
                         ifelse(method == "efron", 1L, 0L))

        ## get results
        deviance <- -2 * (tmp$loglik[1] - tmp$loglik[2])
        covMat <- tmp$imat
        betas <- tmp$coef
        
    } else {
        ## fit the Cox model
        tmp <- survival::coxph(Surv(survTimes, censInd) ~ X, ties=method)

        ## get results
        deviance <- -2 * (tmp$loglik[1] - tmp$loglik[2])
        covMat <- tmp$var
        betas <- tmp$coefficients
    }

    ## return the results
    return(list(betas=betas,
                cov=covMat,
                deviance=deviance))
}
