% Please edit documentation in R/summary.R
\docType{methods}
\name{summary,heart-method}
\alias{summary,heart-method}
\title{Summarizing HEART Model Fits}
\usage{
\S4method{summary}{heart}(object, showcall = TRUE, showpieces = TRUE,
  digits = 3)
}
\arguments{
\item{object}{heart object from \code{heart}.}

\item{showcall}{logic value (TRUE or FALSE) with dafault as TRUE,
indicating whether method \code{show} for object summary.heart prints out
the call information of original call of \code{heart}.}

\item{showpieces}{logic value (TRUE or FALSE) with default as TRUE,
indicating whether method \code{show} for object summary.heart prints out
the baseline pieces.}

\item{digits,}{an interger specifying the desired number of decimal places
(round) for estimates. Negative values are allowed
(\code{help(round)} for more details).}
}
\value{
summary.heart object
}
\description{
\code{summary} returns summary of estimates from HEART model.
}
\details{
These are details.
}
\examples{
data(simuDat)
heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, data = simuDat, 
    baselinepieces = seq(28, 168, length = 5))
summary(heartfit)
}
