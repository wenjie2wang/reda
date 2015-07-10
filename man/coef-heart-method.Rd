% Please edit documentation in R/coef_heart.R
\docType{methods}
\name{coef,heart-method}
\alias{coef,heart-method}
\title{Fitting Heart Model: Piece-wise Gamma Frailty Model for Recurrent Events.
The model is named after the paper title of \emph{Fu et al. (2014)},
Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models}
\usage{
\S4method{coef}{heart}(object, ...)
}
\arguments{
\item{...}{Numeric, complex, or logical vectors.}

\item{na.rm}{A logical scalar.}
}
\value{
If all inputs are integer and logical, then the output
  will be an integer. If integer overflow
  \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
  will be NA with a warning. Otherwise it will be a length-one numeric or
  complex vector.
}
\description{
\code{functionname} returns fitted model results.
}
\details{
This is a test Roxygen comments
}
\examples{
sum(1:10)
sum(1:5, 6:10)
sum(F, F, F, T, T)
function coef for heart object
}
