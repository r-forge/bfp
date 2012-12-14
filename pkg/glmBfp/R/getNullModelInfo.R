#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Time-stamp: <[getNullModelInfo.R] by DSB Mit 21/11/2012 10:00 (CET)>
##
## Description:
## Internal helper function to compute all information about the null model.
##
## History:
## 17/02/2010   file creation
## 15/03/2010   Remove dead comments and also the logGPrior argument from the function.
## 25/05/2010   Remove inclusion of utilities.R
## 28/06/2010   switch back to the Laplace approximation, because it is much
##              more stable for models where there is only a small region where
##              the likelihood is not zero. The "integrate" routine can fail easily
##              there!
## 30/06/2010   include the 2 * pi constant in the log marginal likelihood
## 21/11/2012   include null deviance in return list for TBFs ("nullDeviance")
#####################################################################################


##' Internal helper function to compute all information about the null model
##'
##' This function is used by \code{\link{glmBayesMfp}}.
##' 
##' @param family the family object including the loglikelihood etc.
##' @return a list including the elements
##' \describe{
##' \item{logMargLik}{log marginal likelihood of the null model}
##' \item{nullDeviance}{the deviance of the null model}
##' }
##' 
##' @keywords internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getNullModelInfo <- function(family)
{
    ## this will be the return list:
    ret <- list()

    ## first we need a vectorized loglikelihood function, which will give us
    ## for each of its mu arguments the corresponding loglikelihood, assuming one such mu
    ## for all observations - as is the case in the null model.
    loglikVectorized <- Vectorize(family$loglik)

    ## final likelihood function of the intercept
    loglikIntercept <- function(beta0, log=TRUE)
    {
        ret <- loglikVectorized(family$linkinv(beta0))
        if(! log)
            ret <- exp(ret)

        ret
    }   
    
    ## compute the marginal likelihood of the null model with only the
    ## intercept, that is: the integral wrt to the intercept of the likelihood.    
    ## compute the log marginal likelihood with the Laplace approximation
    ## to be in line with the other marginal likelihood approximations:

    ## optimize
    gaussApproxOptim <- optim(par=mean(family$linPredStart) - 0.1, # add a small
                                        # perturbation deliberately
                              fn=loglikIntercept,
                              method="BFGS",
                              hessian=TRUE,
                              control=
                              list(fnscale=-1))
    
    ## the mode
    beta0Mode <- gaussApproxOptim$par

    ## the precision
    beta0Prec <- - gaussApproxOptim$hess

    ## value of the loglikelihood at the mode
    modeLogLik <- gaussApproxOptim$val

    ## also the null deviance is given by
    ret$nullDeviance <- - 2 * modeLogLik

    ## so the Laplace approximation is
    ret$logMargLik <- modeLogLik + 0.5 * log(2 * pi / beta0Prec)    

    ## --------------------------------------------------
    ## return the list
    return(ret)
}
