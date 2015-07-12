% Please edit documentation in R/coef.R
\docType{methods}
\name{coef,heart-method}
\alias{coef,heart-method}
\title{Extract coefficients estiamtes from HEART model.}
\usage{
coef(object)
}
\arguments{
\item{object}{heart object.}

\item{...}{other arguments.}
}
\value{
A named numeric vector.
}
\description{
\code{coef} is a S4 class generic function
which extracts model coefficients from the heart object returned by modeling
function \code{heart}.
}
\examples{
data(simuDat)
heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, data = simuDat, 
    baselinepieces = seq(28, 168, length = 5))
coef(heartfit)
}
