\name{forecTheta-package}
\alias{forecTheta-package}
\alias{forecTheta}
\docType{package}
\title{
	Forecasting Time Series by Theta Method
	%%~~ package title ~~
}
\description{
	In this package we implemented functions for forecast univariate time series using the Optimised Theta Method of Fioruci et al (2015) and 
	the Standard Theta Method of Assimakopoulos \& Nikolopoulos (2000). 
	Moreover, it is including a function for compute the main errors metrics used in time series forecasting and a function for compute the 
	Generalised Rolling Origin Evaluation, which contain as particular cases the Rolling Origin Evaluation and the Fixed Origin Evaluation of Tashman (2000).
	%%~~ A concise (1-5 lines) description of the package ~~
}
\details{
	\tabular{ll}{
		Package: \tab forecTheta\cr
		Type: \tab Package\cr
		Version: \tab 1.0\cr
		Date: \tab 2015-03-30\cr
		License: \tab GPL (>=2.0) \cr
	}
	
	%%~~ An overview of how to use the package, including the most important ~~
	%%~~ functions ~~
	\code{otm(y, h, seasonal=FALSE, theta=NULL)}
	
	\code{thetaM(y, h, seasonal=FALSE)}
	
	\code{errorMetric(obs, forec, type = "sAPE", statistic = "M")}
	
	\code{groe(y, forecFunction = ses, g = "sAPE", n1 = length(y)-10)}
}
\author{
	Jose Augusto Fioruci, Francisco Louzada and Bao Yiqi

	Maintainer: Jose Augusto Fioruci <jafioruci@gmail.com>
	%%~~ The author and/or maintainer of the package ~~
}
\references{
	Fioruci J.A., Pellegrini T.R., Louzada F., Petropoulos F. (2015). \emph{The Optimised Theta Method}. Free available at \url{ http://arxiv.org/abs/1503.03529 }. 
	
	Assimakopoulos, V. and Nikolopoulos k. (2000). \emph{The theta model: a decomposition approach to forecasting}. International Journal of Forecasting 16, 4, 521--530.

	Tashman, L.J. (2000). \emph{Out-of-sample tests of forecasting accuracy: an analysis and review}. International Journal of Forecasting 16 (4), 437--450.
	%%~~ Literature or other references for background information ~~
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%%~~ the R documentation directory ~~
\keyword{ Optimised Theta Method, Theta Method, Rolling Origin Evaluation, Fixed Origin Evaluation, Time Series Forecasting, Error Metric,  Cross Validation  }
\seealso{
	\code{\link[forecTheta]{otm}}, \code{\link[forecTheta]{thetaM}}, \code{\link[forecTheta]{groe}}, 
	\code{\link[forecTheta]{rolOrig}}, \code{\link[forecTheta]{fixOrig}},
	\code{\link[forecTheta]{errorMetric}}
	%%~~ Optional links to other man pages, e.g. ~~
	%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%% ~~ simple examples of the most important functions ~~
##############################################################	

y1 = 2+ 0.15*(1:20) + rnorm(20,2)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30,2)
y =  as.ts(c(y1,y2))

otm.fit <- otm(y=as.ts(y[1:40]), h=10)
plot(y)
points(otm.fit$fitted, type = "l")
points(otm.fit$mean, type = "l", col="blue")

theta.fit <- thetaM(y=as.ts(y[1:40]), h=10)

### sMAPE metric
errorMetric(obs=as.ts(y[41:50]), forec=otm.fit$mean, type = "sAPE", statistic = "M")
errorMetric(obs=as.ts(y[41:50]), forec=theta.fit$mean, type = "sAPE", statistic = "M")

### sMdAPE metric
errorMetric(obs=as.ts(y[41:50]), forec=otm.fit$mean, type = "sAPE", statistic = "Md")
errorMetric(obs=as.ts(y[41:50]), forec=theta.fit$mean, type = "sAPE", statistic = "Md")

### MASE metric
meanDiff1 = mean(abs(diff(as.ts(y[1:40]), lag = 1)))
errorMetric(obs=as.ts(y[41:50]), forec=otm.fit$mean, type = "AE", statistic = "M") / meanDiff1
errorMetric(obs=as.ts(y[41:50]), forec=theta.fit$mean, type = "AE", statistic = "M") / meanDiff1
}
