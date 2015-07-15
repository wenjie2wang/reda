% Please edit documentation in R/plot.R
\name{plot_sampleMCF}
\alias{plot_sampleMCF}
\title{Plot Empirical Mean Cumulative Function (MCF) from Sample Data.}
\usage{
plot_sampleMCF(MCF, linetypes, linecolors, ...)
}
\arguments{
\item{MCF}{data.frame generated from function
\code{\link[heart]{sample_MCF}}.}

\item{linetypes}{line types specified to different groups.}

\item{linecolors}{line colors specified to different groups.}

\item{...}{further arguments.}
}
\value{
ggplot object.
}
\description{
\code{plot_sampleMCF} plots empirical mean cumulative function (MCF).
}
\details{
The function plots empirical mean cumulative function
by using ggplot2 plotting system.
So the plots generated are able to be further customized properly.
}
