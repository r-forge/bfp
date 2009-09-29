\name{ozone}
\alias{ozone}
\docType{data}
\title{Ozone data from Breiman and Friedman, 1985}
\description{
This is the Ozone data discussed in Breiman and Friedman (JASA, 1985, p. 580).
These data are for 330 days in 1976.  All measurements are in the area of
Upland, CA, east of Los Angeles.
}
\format{
  A data frame with 366 observations on the following 13 variables.
  \describe{
    \item{\code{month}}{a numeric vector}
    \item{\code{day}}{a numeric vector}
    \item{\code{weekday}}{a factor with levels \code{Monday} \code{Tueday} \code{Wednesday} \code{Thursday} \code{Friday} \code{Saturday} \code{Sunday}}
    \item{\code{hourAverageMax}}{a numeric}
    \item{\code{pressure500Height}}{a numeric}
    \item{\code{windSpeed}}{a numeric}
    \item{\code{humidity}}{a numeric}
    \item{\code{tempSandburg}}{a numeric}
    \item{\code{tempElMonte}}{a numeric}
    \item{\code{inversionBaseHeight}}{a numeric}
    \item{\code{pressureGradientDaggett}}{a numeric}
    \item{\code{inversionBaseTemp}}{a numeric}
    \item{\code{visibility}}{a numeric}
  }
}


\usage{data(ozone)}
\source{Breiman, L and  Friedman, J.  (1985),
 \dQuote{Estimating Optimal Transformations for Multiple 
Regression and Correlation},
 \emph{Journal of the American Statistical Association}, 80,  580-598.
    }
\keyword{datasets}
