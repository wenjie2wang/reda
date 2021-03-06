% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R
\docType{class}
\name{rateReg-class}
\alias{rateReg-class}
\title{An S4 Class Representing a Fitted Model}
\description{
The class \code{rateReg} is an S4 class that represents a fitted model.  The
function \code{\link{rateReg}} produces objects of this class.  See
``Slots'' for details.
}
\section{Slots}{

\describe{
\item{\code{call}}{Function call.}

\item{\code{formula}}{Formula.}

\item{\code{nObs}}{A positive integer}

\item{\code{spline}}{A list.}

\item{\code{estimates}}{A list.}

\item{\code{control}}{A list.}

\item{\code{start}}{A list.}

\item{\code{na.action}}{A character vector (of length one).}

\item{\code{xlevels}}{A list.}

\item{\code{contrasts}}{A list.}

\item{\code{convergCode}}{A nonnegative integer.}

\item{\code{logL}}{A numeric value.}

\item{\code{fisher}}{A numeric matrix.}
}}

\seealso{
\code{\link{rateReg}}
}
