#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Course:
##         from the Department of Statistics, University of Munich
## Time-stamp: <[test.R] by DSB Fre 04/12/2009 17:38 (CET)>
##
## Description:
##
##
## History:
## 29/11/2008   update for new model prior option
## 22/09/2009   update for new SWITCH move type: use 2 FP terms in order to be able
##              to do the simulation
#####################################################################################

library(bfp)

## setting

beta0 <- 1
alpha1 <- 1
alpha2 <- 3
delta1 <- 1

sigma <- 2                              # sigma2 = 4
n <- 15
k <- 2L

## simulate data

set.seed (123)

x <- matrix (runif (n * k, 1, 4), nrow = n, ncol = k) # predictor values
w <- matrix (rbinom (n * 1, size = 1, prob = 0.5), nrow = n, ncol = 1)

x1tr <- alpha1 * x[,1]^2
x2tr <- alpha2 * (x[,2])^(1/2)
w1tr <- delta1 * w[,1]

predictorTerms <-
    x1tr +
    x2tr +
    w1tr

trueModel <- list (powers = list (x1 = 2, x2 = 0.5),
                   ucTerms = as.integer (1)
                   )

covariateData <- data.frame (x1 = x[,1],
                             x2 = x[,2],
                             w = w)

covariateData$y <- predictorTerms + rnorm (n, 0, sigma)
covariateData

exhaustive <- BayesMfp (y ~ bfp (x1, max=1) + bfp(x2, max=1),
                        data = covariateData,
                        priorSpecs =
                        list (a = 3.5,
                              modelPrior="flat"),
                        method = "exhaustive",
                        nModels = 100
                        )           
attr(exhaustive, "logNormConst")
summary(exhaustive)

truedata <- as.data.frame(exhaustive)
truedata

post <- exp(truedata$logMargLik + truedata$logPrior)
normConst <- sum(post)
post / normConst

set.seed(2)
simulation <- BayesMfp (y ~ bfp (x1, max=1) + bfp(x2, max=1),
                        data = covariateData,
                        priorSpecs =
                        list (a = 3.5,
                              modelPrior="flat"),
                        method = "sampling",
                        nModels = 100,
                        chainlength=1000000
                        )           
attr(simulation, "logNormConst")
summary(simulation)

p1 <- posteriors(exhaustive)
p2 <- posteriors(simulation, 2)

plot(log(p1), log(p2))
abline(0, 1)

d1 <- as.data.frame(exhaustive)
d2 <- as.data.frame(simulation)

head(d1[, -c(7,8)])
head(d2[, -c(2, 8, 9)])

shouldbeZero <- d1[, -c(7,8)] - d2[, -c(2, 8, 9)]
shouldbeZero <- max(abs(unlist(shouldbeZero)))
stopifnot(all.equal(shouldbeZero, 0))

