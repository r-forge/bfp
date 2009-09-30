#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[rshrinkage.R] by DSB Fre 05/09/2008 13:17 (CEST) on daniel@puc.home>
##
## Description:
## Sample from the model-specific posterior of the shrinkage factor t = g / (1 + g).
##
## History:
## 04/09/2008   file creation
## 05/09/2008   correct innerFraction computation
#####################################################################################


rshrinkage <- function(n=1,             # number of samples
                       R2,              # coefficient of determination in the model
                       nObs,            # number of observations
                       p,               # number of effects (without intercept)
                       alpha            # hyperparameter for hyper-g prior
                       )
{
    ## uniforms needed for inverse sampling
    uniforms <- runif(n=n)

    ## parameters for the beta functions
    shape1 <- (nObs - p - alpha + 1) / 2
    shape2 <- (p + alpha - 2) / 2
    
    ## evaluate the quantile function at these values
    oneMinusR2 <- 1 - R2

    innerFraction <- (1 - R2) / qbeta(uniforms + (1 - uniforms) * pbeta(1 - R2, shape1, shape2),
                                      shape1, shape2)
              
    ret <- (1 - innerFraction) / R2

    return(ret)    
}
