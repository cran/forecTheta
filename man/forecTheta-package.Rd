\name{forecTheta-Package}
\alias{forecTheta-package}
\alias{forecTheta}
\docType{package}
\title{
	Forecasting Time Series by Theta Models
	%%~~ package title ~~
}
\description{
This package provides functions for forecasting univariate time series using several Theta models, originally proposed by Assimakopoulos and Nikolopoulos (2000) and later extended by Fiorucci et al. (2016). This version also includes implementations of bagging methods, based on the work of Bergmeir et al. (2016), applied to the DOTM, DSTM, OTM, and STM models.
	%%~~ A concise (1-5 lines) description of the package ~~
}
\details{
	\tabular{ll}{
		Package: \tab forecTheta\cr
		Type: \tab Package\cr
		Version: \tab 3.0\cr
		Date: \tab 2025-05-20\cr
		License: \tab GPL (>=2.0) \cr
	}

	%%~~ An overview of how to use the package, including the most important ~~
	%%~~ functions ~~
	\code{dotm(y, h)}

	\code{bagged_dotm(y, h)}

	\code{stheta(y, h)}

	\code{errorMetric(obs, forec, type = "sAPE", statistic = "M")}

	\code{PI_eval(obs, forec, lower_bounds, upper_bounds, name = c("MSIS", "ACD"), alpha = 0.05)}

	\code{groe(y, forecFunction = ses, g = "sAPE", n1 = length(y)-10)}
}
\author{
	Jose Augusto Fiorucci, Francisco Louzada, Igor De Oliveira Barros Faluhelyi.

	Maintainer: Jose Augusto Fiorucci <jafiorucci@gmail.com>
	%%~~ The author and/or maintainer of the package ~~
}
\references{
  Assimakopoulos, V. and Nikolopoulos k. (2000). \emph{The theta model: a decomposition approach to forecasting}. International Journal of Forecasting 16, 4, 521--530, <doi:10.1016/S0169-2070(00)00066-2>.

	Bergmeir, C., Hyndman, R.J. and Benítez, J. M. (2016). \emph{Bagging exponential smoothing methods using STL decomposition and Box–Cox transformation}. International journal of forecasting 32 (2), 303--312, <doi:10.1016/j.ijforecast.2015.07.002>.

	Fiorucci J.A., Pellegrini T.R., Louzada F., Petropoulos F., Koehler, A. (2016). \emph{Models for optimising the theta method and their relationship to state space models},	International Journal of Forecasting, 32 (4), 1151--1161, <doi:10.1016/j.ijforecast.2016.02.005>.

	Fioruci J.A., Pellegrini T.R., Louzada F., Petropoulos F. (2015). \emph{The Optimised Theta Method}. arXiv preprint, arXiv:1503.03529.

	Gneiting, T. and Raftery, A.E. (2007). Strictly proper scoring rules, prediction, and estimation. Journal of the American Statistical Association 102 (477), 359--378, \doi{10.1198/016214506000001437}.

	Tashman, L.J. (2000). \emph{Out-of-sample tests of forecasting accuracy: an analysis and review}. International Journal of Forecasting, 16 (4), 437--450, <doi:10.1016/S0169-2070(00)00065-0>.
	%%~~ Literature or other references for background information ~~
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%%~~ the R documentation directory ~~
\keyword{Dynamic Optimised Theta Model}
\keyword{Theta Method}
\keyword{time series forecasting}




\seealso{
	\code{\link[forecTheta]{dotm}}, \code{\link[forecTheta]{bagged_dotm}}, \code{\link[forecTheta]{stheta}}, \code{\link[forecTheta]{otm.arxiv}}, \code{\link[forecTheta]{groe}}, \code{\link[forecTheta]{rolOrig}}, \code{\link[forecTheta]{fixOrig}},
	\code{\link[forecTheta]{errorMetric}}
	%%~~ Optional links to other man pages, e.g. ~~
	%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%% ~~ simple examples of the most important functions ~~
##############################################################

y1 = 2+ 0.15*(1:20) + rnorm(20)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30)
y =  as.ts(c(y1,y2))
out <- dotm(y, h=10)
summary(out)
plot(out)

out <- dotm(y=as.ts(y[1:40]), h=10)
summary(out)
plot(out)

out2 <- bagged_dotm(y=as.ts(y[1:40]), h=10)
summary(out2)
plot(out2)

out3 <- stheta(y=as.ts(y[1:40]), h=10)
summary(out3)
plot(out3)

### sMAPE metric
errorMetric(obs=as.ts(y[41:50]), forec=out$mean, type = "sAPE", statistic = "M")
errorMetric(obs=as.ts(y[41:50]), forec=out2$mean, type = "sAPE", statistic = "M")
errorMetric(obs=as.ts(y[41:50]), forec=out3$mean, type = "sAPE", statistic = "M")

### sMdAPE metric
errorMetric(obs=as.ts(y[41:50]), forec=out$mean, type = "sAPE", statistic = "Md")
errorMetric(obs=as.ts(y[41:50]), forec=out2$mean, type = "sAPE", statistic = "Md")
errorMetric(obs=as.ts(y[41:50]), forec=out3$mean, type = "sAPE", statistic = "Md")

### MASE metric
meanDiff1 = mean(abs(diff(as.ts(y[1:40]), lag = 1)))
errorMetric(obs=as.ts(y[41:50]), forec=out$mean, type = "AE", statistic = "M") / meanDiff1
errorMetric(obs=as.ts(y[41:50]), forec=out2$mean, type = "AE", statistic = "M") / meanDiff1
errorMetric(obs=as.ts(y[41:50]), forec=out3$mean, type = "AE", statistic = "M") / meanDiff1

}
