#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Course:
##         from the Department of Statistics, University of Munich
## Time-stamp: <[testSimulation.R] by DSB Sam 29/11/2008 19:12 (CET) on daniel@puc.home>
##
## Description:
##
##
## History:
## 14/10/2008   file creation to detect the error leading to an incompatible
##              dimensions exception
## 29/11/2008   update for new model prior option
#####################################################################################

library(bfp)

set.seed (234)

## setting where the error occured
beta0 <- 1
alpha1 <- 1
alpha2 <- c (1, 1)
delta1 <- 1

sampleSize <- c (40, 100)
sigma2 <- c (4, 3, 2, 1)
hyperpara <- c (3.05, 3.4, 3.7, 3.95)

h <- 1
i <- 4
j <- 4

thisN <- sampleSize[h]

x <- matrix (runif (thisN * 3, 1, 4), nrow = thisN, ncol = 3) # predictor values
w <- matrix (rbinom (thisN * 2, size = 1, prob = 0.5), nrow = thisN, ncol = 2)

covData <- data.frame (x = x, w = w)# start data frame

x1tr <- alpha1 * x[,1]^2
x2tr <- cbind (x[,2]^-0.5, x[,2]^-0.5 * log (x[,2])) %*% alpha2
w1tr <- delta1 * w[,1]

predictorTerms <- x1tr + x2tr + w1tr # linear predictor
    
thisPriorSpecs <- list(a = hyperpara[i],
                       modelPrior="sparse")
        
covData$y <- predictorTerms + rnorm (thisN, 0, sqrt (sigma2[j]))

## try it:
modelNow <- BayesMfp (y ~ bfp (x.1) + bfp (x.2) + bfp (x.3) + uc (w.1) + uc (w.2),
                      data = covData,
                      priorSpecs = thisPriorSpecs,
                      method = "exhaustive"
                      )
