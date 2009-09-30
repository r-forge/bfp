#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[simulation.R] by DSB Sam 29/11/2008 19:03 (CET) on daniel@puc.home>
##
## Description:
## Example for exhaustive model enumeration.
##
## History:
## 04/07/2008   copy from full simulation example.
## 29/11/2008   use default values for priorSpecs
#####################################################################################

## small simulation study

## setting

beta0 <- 1
alpha1 <- 1
alpha2 <- c (1, 1)
delta1 <- 1

sigma <- 2                              # sigma2 = 4
n <- 40

## simulate data

set.seed (123)

x <- matrix (runif (n * 3, 1, 4), nrow = n, ncol = 3) # predictor values
w <- matrix (rbinom (n * 2, size = 1, prob = 0.5), nrow = n, ncol = 2)

x1tr <- alpha1 * x[,1]^2
x2tr <- cbind (x[,2]^-0.5, x[,2]^-0.5 * log (x[,2])) %*% alpha2

summary (x1tr / x2tr)
mean (x1tr) / mean (x2tr)

w1tr <- delta1 * w[,1]

predictorTerms <-
    x1tr +
    x2tr +
    w1tr

trueModel <- list (powers = list (x.1 = 2, x.2 = c (-0.5, -0.5), x.3 = numeric (0)),
                   ucTerms = as.integer (1)
                   )

covariateData <- data.frame (x = x, w = w)

## one simulation:

covariateData$y <- predictorTerms + rnorm (n, 0, sigma)

system.time (
             simulation1Diagonal <- BayesMfp (y ~ bfp (x.1) + bfp (x.2) + bfp (x.3) + uc (w.1) + uc (w.2),
                                              data = covariateData,
                                              method = "exhaustive"
                                              )
             )
ind <- findModel (trueModel, simulation1Diagonal)

estimate <- simulation1Diagonal[ind]
summary (estimate)

plotx1 <- plotCurveEstimate (estimate, "x.1")
with (plotx1, lines (original, original^2, col = "red"))
plotx2 <- plotCurveEstimate (estimate, "x.2")
with (plotx2, lines (original, original^-0.5 * (1 + log (original)), col = "red"))

simul1Data <- as.data.frame (simulation1Diagonal)
simul1Data$rankML <- as.character (rank (-simul1Data$logMargLik))
simul1Data[1:20, ]

plotCurveEstimate (simulation1Diagonal[1], "x.1")

## save (simulation1Diagonal, simul1Data, file = c ("simulation1Diagonal2.RData"))
## load ("simulation1Diagonal2.RData")

library (doBy)
orderBy (~ -logMargLik, data = simul1Data)[1:100,] # order by marginal likelihood
