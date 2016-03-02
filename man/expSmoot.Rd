\name{expSmoot}
\alias{expSmoot}

\title{
	Simple Exponential Smoothing Method
}
\description{
	Estimation of Simple Exponential Smoothing Method
}
\usage{
	expSmoot(y, h=5, ell0=NULL, alpha=NULL, lower = c(-1e+10, 0.1), 
	upper = c(1e+10, 0.99))
}

\arguments{
	\item{y}{
		Object of time series class.
	}
	\item{h}{
		Number of required forecasting periods.
	}
	\item{ell0}{
		The value of \code{ell0^*} parameter.
	}
	\item{alpha}{
		The value of \code{alpha} parameter.
	}
	\item{lower}{
		The lower limit of parametric space.
	}
	\item{upper}{
		The upper limit of parametric space.
	}
}

\value{
	A list containing the elements:
	\item{$y }{The original time series.}
    \item{$par }{The estimated values for \code{(ell^*, alpha)} parameters}
	\item{$mean}{The forecasting values}
    \item{$fitted }{A time series element with the fitted points.}
    \item{$residuals }{A time series element with the residual points.}
}

\author{
	Jose Augusto Fiorucci, Francisco Louzada and Bao Yiqi
}

\seealso{
	\code{\link[forecTheta]{forecTheta-package}}, \code{\link[forecTheta]{otm.arxiv}}, \code{\link[forecTheta]{dotm}}
}

\examples{

y1 = 2+ 0.15*(1:20) + rnorm(20,2)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30,2)
y =  as.ts(c(y1,y2))

expSmoot(y, h=10)

}

\keyword{simple exponential smoothing}

