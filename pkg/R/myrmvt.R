#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[myrmvt.R] by DSB Fre 04/07/2008 11:37 (CEST) on daniel@puc.home>
##
## Description:
## Helper function for sampling from multivariate normal distribution.
##
## History:
## 04/07/2008   copy from thesis function collection.
#####################################################################################

myrmvt <- function (n, sigma = diag (2), mu = rep (0, 2), df = 1)
{
    ret <- rmvt (n = n, sigma = sigma, df = df)
    ret <- sweep (ret, 2, mu, "+")
    ret
}
