% Please edit documentation in R/coef.R
\docType{methods}
\name{confint,heart-method}
\alias{confint,heart-method}
\title{Confidence Intervals for HEART Model Coefficients}
\usage{
confint(object, parm, level = 0.95, ...)
}
\arguments{
\item{object}{heart object.}

\item{parm}{a specification of which parameters are
to be given confidence intervals,
either a vector of numbers or a vector of names.
If missing, all parameters are considered.}

\item{level}{the confidence level required.}

\item{...}{additional argument(s).}
}
\value{
A numeric matrix with rownames and colnames.
}
\description{
\code{confint} is a S4 class generic function for heart object,
which computes confidence intervals
for all or selected coefficients in a fitted HEART model.
}
\details{
Under regularity condition (Shao, 2003,
Theorem 4.16 and Theorem 4.17, page 287, 290),
the confidence intervals are constructed loosely
based on Fisher information matrix and estimates of coefficients.
See \emph{Fu et al. (2014)} for more details.
}
\examples{
data(simuDat)
heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, data = simuDat, 
    baselinepieces = seq(28, 168, length = 5))
confint(heartfit)
}
\references{
Shao, J. (2003),
\emph{Mathematical statistics}, Springer texts in statistics,
New York: Springer, 2nd edition.

Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
"Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models,"
\emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
}
