% Please edit documentation in R/MCF.R
\name{sample_MCF}
\alias{sample_MCF}
\title{Empirical Mean Cumulative Function (MCF)}
\usage{
sample_MCF(formula, data, subset, na.action, ...)
}
\arguments{
\item{formula}{Survr object from function \code{Survr} in package survrec.
The covariate specified in the rhs of the formula can either be 1 or
any one factor variable in the data.  The former computes the overall
empirical MCF from sample.  The latter computes the empirical MCF for each
level of the factor variable specified respectively.}

\item{data}{an optional data frame, list or environment containing
the variables in the model.  If not found in data, the variables are taken
from \code{environment(formula)}, usually the environment from which
the function is called.}

\item{subset}{an optional vector specifying a subset of observations
to be used in the fitting process.}

\item{na.action}{a function which indicates what should the procedure do
if the data contains NAs.  The default is set by the
na.action setting of \code{\link[base]{options}} and is na.fail if that is
not set.  The "factory-fresh" default is \code{\link[stats]{na.omit}}.
Another possible value is NULL, no action.
Value \code{\link[stats]{na.exclude}} can be useful.}

\item{...}{further arguments}
}
\value{
data.frame
}
\description{
\code{sample_MCF} computes empirical mean cumulative function (MCF) on every
time point from recurrent event sample data.
It does not assume any particular underlying model.
}
