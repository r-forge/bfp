


\name{SMARTfull}
\alias{SMARTfull}
\docType{data}
\title{Imputed and transformed data extract from the SMART study}
\description{
Data from the SMART study, a prospective cohort study at the University
Medical Centre Utrecht, the Netherlands. The study was designed to (a) establish the
prevalence of concomitant arterial diseases and risk factors for cardiovascular 
disease in a high-risk population; (b) identify predictors of future cardiovascular 
events in patients with symptomatic cardiovascular disease.
}
\usage{data("SMARTfull")}
\format{
  A data frame with 3873 observations on the following 16 variables.
  \describe{
    \item{\code{TEVENT}}{a numeric vector. Time to cardiovascular event (days)}
    \item{\code{EVENT}}{a logical vector. Cardiovascular event (clinical, 0/1)}
    \item{\code{AGE.TRANS}}{a numeric vector. (Age(years)-50)^2 if over 50, otherwise 0}
    \item{\code{SEX}}{a numeric vector. \code{1} Male \code{2} Female}
    \item{\code{SMOKING}}{a factor with levels \code{1} Never \code{2} Former \code{3} Current}
    \item{\code{ALCOHOL}}{a factor with levels \code{1} Never \code{2} Former \code{3} Current}
    \item{\code{BMI}}{a numeric vector (imputed). Body mass index (kg/m^2)}
    \item{\code{SYSTH}}{a numeric vector (imputed).  Systomlic blood pressure (mm/Hg)}
    \item{\code{HDL}}{a numeric vector (imputed). High-density lipoprotein cholesterol (mmol/L)}
    \item{\code{DIABETES}}{a numeric vector (imputed). Ever diabetes (0/1)}
    \item{\code{HISTCAR2}}{a numeric vector. History of cardiac disease score.}
    \item{\code{HOMOC}}{a numeric vector (imputed). Homocysteine (\eqn{\mu}mol/L}
    \item{\code{logCREAT}}{a numeric vector (imputed). Creatinine clearance (log mL/min)}
    \item{\code{ALBUMIN}}{a factor with levels \code{1} No \code{2} Micro \code{3} Macro. Albumin in urine}
    \item{\code{STENOSIS}}{a numeric vector (imputed). Carotic stenosis \eqn{\ge 50\%} by duplex}
    \item{\code{IMT}}{a numeric vector (imputed). Intima media thickness (mm)}
  }
}
\details{
\code{SMARTfull} is a processed extract from the SMART study. The code for processing is available: \url{http://www.clinicalpredictionmodels.org/doku.php?id=rcode_and_data:chapter23}. This dataset is named \code{SMARTc} in that code, and is used to fit the "full" model, which is discussed in Chapter 23 of Steyerberg's \emph{Clinical Prediction Models}.}

\source{\url{http://www.clinicalpredictionmodels.org}}
\references{Steyerberg, E.W. (2009) \emph{Clinical Prediction Models}. New York: Springer-Verlag }
\examples{
data(SMARTfull)
}
\keyword{datasets}
