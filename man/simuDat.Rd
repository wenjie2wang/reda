% Please edit documentation in R/data.R
\docType{data}
\name{simuDat}
\alias{simuDat}
\title{Simulated Dataset for Demonstration}
\format{data frame}
\description{
A data frame with covariates named
'ID', 'time', 'event', 'group' and 'X1'.
}
\details{
The dataset is simulated by the thinning method developed
by \emph{Lewis and Shedler (1979)}.
See \emph{Fu et al. (2014)} for more details.
}
\examples{
data(simuDat)
str(simuDat)
}
\references{
Lewis, P. and Shedler, G. (1979),
“Simulation of nonhomogeneous Poisson processes by thinning,”
\emph{Naval Research Logistics Quarterly}, 26, 403–413.

Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
"Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models,"
\emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
}
