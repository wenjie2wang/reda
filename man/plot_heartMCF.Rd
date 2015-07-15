% Please edit documentation in R/plot.R
\name{plot_heartMCF}
\alias{plot_heartMCF}
\title{Plot Estimated Mean Cumulative Function (MCF) of Baseline Rate Function.}
\usage{
plot_MCF(object, CI = TRUE, level = 0.95, ...)
}
\arguments{
\item{object}{heart object.}

\item{CI}{logic value, TRUE or FALSE to specify whether to include confidence
interval in the plot.}

\item{level}{numeric value \eqn{\in (0, 1)} confidence level.}

\item{...}{further arguments.}
}
\value{
ggplot object
}
\description{
\code{plot_heartMCF} produces the plot of baseline function.
}
\details{
The function estimates mean cumulative function of baseline rate function
by using ggplot2 plotting system.
So the plots generated are able to be further customized properly.
}
