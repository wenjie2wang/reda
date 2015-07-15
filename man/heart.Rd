% Please edit documentation in R/fit.R, R/heart.R
\docType{package}
\name{heart}
\alias{heart}
\alias{heart-package}
\title{Fitting HEART Model}
\usage{
heart(formula, baselinepieces, data, subset, na.action, start = list(), 
    control = list(), contrasts = NULL, ...)
}
\arguments{
\item{formula}{Survr object from function \code{Survr} in package survrec.}

\item{baselinepieces}{an optional numeric vector consisting of
all the right endpoints of baseline pieces.}

\item{data}{an optional data frame, list or environment containing
the variables in the model.  If not found in data, the variables are taken
from \code{environment(formula)}, usually the environment from which
\code{heart} is called.}

\item{subset}{an optional vector specifying a subset of observations
to be used in the fitting process.}

\item{na.action}{a function which indicates what should the procedure do
if the data contains NAs.  The default is set by the
na.action setting of \code{\link[base]{options}} and is na.fail if that is
not set.  The "factory-fresh" default is \code{\link[stats]{na.omit}}.
Another possible value is NULL, no action.
Value \code{\link[stats]{na.exclude}} can be useful.}

\item{start}{an optional list of starting values for the parameters
to be estimated in the model.}

\item{control}{an optional list of parameters for controlling the likelihood
function maximization process. For more detail, users may \code{help(nlm)}.}

\item{contrasts}{an optional list, whose entries are values
(numeric matrices or character strings naming functions) to be used
as replacement values for the contrasts replacement function and
whose names are the names of columns of data containing factors.
See the \code{contrasts.arg} of \code{model.matrix.default} for more detail.}
}
\value{
a heart object.
}
\description{
\code{heart} returns fitted model results.

The heart package provides functions to fit HEART model for recurrent event
analysis.
}
\details{
HEART model is a piece-wise Gamma frailty model for recurrent events.
The model is named after the paper title of \emph{Fu et al. (2014)},
Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models
}
\examples{
data(simuDat)
heartfit <- heart(formula = Survr(ID, time, event) ~ X1 + group, data = simuDat, 
    baselinepieces = seq(28, 168, length = 6))
str(heartfit)
show(heartfit)  # or simply call heartfit
summary(heartfit)
coef(heartfit)
confint(heartfit)
baseline(heartfit)
plot_MCF(heartfit)
}
\references{
Fu, Haoda, Junxiang Luo, and Yongming Qu. (2014),
"Hypoglycemic Events Analysis via Recurrent Time-to-Event (HEART) Models,"
\emph{Journal of biopharmaceutical statistics}, 2014 Dec 1, Epub 2014 Dec 1.
}
