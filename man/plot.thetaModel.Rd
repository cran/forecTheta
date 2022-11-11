\name{Plot}
\alias{plot.thetaModel}
\title{Plot forecasts points and prediction intervals for \code{thetaModel} objects}
\usage{\method{plot}{thetaModel}(x, ylim=NULL, xlim=NULL, ylab=NULL, xlab=NULL, main=NULL, ...)
}

\arguments{
\item{x}{Object of class \dQuote{thetaModel}.}
\item{ylim}{the y limits of the plot.}
\item{xlim}{the x limits of the plot.}
\item{ylab}{a label for the y axis.}
\item{xlab}{a label for the x axis.}
\item{main}{a main title for the plot.}
\item{...}{Other plotting parameters passed to \code{\link[graphics]{par}}.}
}

\description{Produces a figure of the time series and the forecasts points from Optimised Theta Method.}

\value{None. Function produces a plot}

\seealso{\code{\link{dotm}}, \code{\link[forecTheta]{forecTheta-package}}}

\author{Jose A Fiorucci}
\examples{
y1 = 2+ 0.15*(1:20) + rnorm(20,2)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30,2)
y =  as.ts(c(y1,y2))
out <- dotm(y, h=10)
plot(out)
}
\keyword{plot}
