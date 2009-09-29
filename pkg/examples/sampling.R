#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[sampling.R] by DSB Fre 04/07/2008 18:14 (CEST) on daniel@puc.home>
##
## Description:
## Example for sampling.
##
## History:
## 04/07/2008   copy from full sampling ozone example.
#####################################################################################

data(ozone)

symnum (cor (subset (na.omit (ozone), select = - weekday))) # correlations between two variables
ozoneFull <- na.omit (subset (ozone, select = - tempElMonte))

workday <- ozoneFull$weekday
levels (workday)
levels (workday) <- c (rep ("workday", 5), rep ("weekend", 2))
ozoneFull$workday <- workday

## day of the year

library (doBy)

ozoneFull <- orderBy (~ month + day, ozoneFull)
ozoneFull$dayOfYear <- seq_len (nf <- nrow (ozoneFull))

set.seed (321)                          # determine test and training set
whichTest <- sample (nf, 30, replace = FALSE)
ozoneTest <- ozoneFull [whichTest, ]
ozoneTraining <- ozoneFull[-whichTest, ]

ozoneTest <- orderBy (~month + day, ozoneTest)
testDates <- with (ozoneTest, paste (month, day, sep = "/"))
testDates <- matrix (testDates, nrow = 3, ncol = 10, byrow = TRUE)


## run sampler
set.seed (432)
system.time (
             ozoneModels <- BayesMfp (hourAverageMax ~
                                      bfp (dayOfYear) +
                                      bfp (pressure500Height) +
                                      bfp (windSpeed) +
                                      bfp (humidity) +
                                      bfp (tempSandburg) +
                                      bfp (inversionBaseHeight) +
                                      bfp (pressureGradientDaggett) +
                                      bfp (inversionBaseTemp) +
                                      bfp (visibility),
                                      data = ozoneTraining,
                                      nModels = 300,
                                      method = "sampling",
                                      chainlength = 300
                                      )
             )

sumOzoneModels <- summary (ozoneModels)
str (sumOzoneModels)
