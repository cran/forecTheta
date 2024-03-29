\name{otm.arxiv}
\alias{otm.arxiv}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Optimised Theta Method
	%%  ~~function to do ... ~~
}
\description{
	%%  ~~ A concise (1-5 lines) description of what the function does. ~~
	Functions for forecast univariate time series using the Optimised Theta Method presented in the arxiv paper (Fioruci et al, 2015).
	If the theta parameter is not specified so the Generalised Rolling Origin Evaluation is used
	for select the theta value over the thetaList argument.
}
\usage{
	otm.arxiv( y, h=5, s=NULL, theta=NULL, tLineExtrap=expSmoot, g="sAPE",
		approach="c", n1=NULL, m=NULL, H=NULL, p=NULL,
		thetaList=seq(from=1,to=5,by=0.5), mc.cores=1, ...)
}

\arguments{
	\item{y}{
		Object of time series class
		%%     ~~Describe \code{y} here~~
	}
	\item{h}{
		Number of required forecasting periods
		%%     ~~Describe \code{h} here~~
	}
	\item{s}{
		If \code{TRUE}, the multiplicative seasonal decomposition is used.
		If \code{NULL}, quarterly and monthly time series are tested for statistically seasonal behaviour, with 95\% of significance. Default is NULL.
		%%     ~~Describe \code{seasonal} here~~
	}
	\item{theta}{
		The value of theta parameter. If \code{theta = NULL} the theta parameter is estimated using the Generalised Rolling Origin Evaluation.
		%% Otherwise, the Optimised Theta Model is computed for \code{theta}.
		%%     ~~Describe \code{theta} here~~
	}
	\item{tLineExtrap}{
		A forecasting function for extrapolation the second theta-line. Default is \code{expSmoot}.
		%%     ~~Describe \code{tLineExtrap} here~~
	}
	\item{g}{
		The error type that will be used by \code{groe} function for select the \code{theta} value in the estimation process.
		The possibility values for \code{g}	is \code{"sAPE", "APE", "AE"} and \code{"SE"}.
		If \code{theta} is not \code{NULL} the \code{g} argument is not used. Default is "sAPE".
		%%     ~~Describe \code{g} here~~
	}
	\item{approach}{
		The approach set-up for \code{groe} parameters (\code{n1, m, H, p}). One letter between \code{'a'} to \code{'h'} according to Fioruci et al (2015).
	}
	\item{n1}{
		The first origin for Generalised Rolling Origin Evaluation.
		This argument is not used if \code{theta!=NULL} or \code{approach!=NULL}.
		%%     ~~Describe \code{n1} here~~
	}
	\item{m}{
		The number of movements of the origin in each step. This argument is not used if \code{theta!=NULL} or \code{approach!=NULL}.
		%%     ~~Describe \code{m} here~~
	}
	\item{H}{
		 The number of predictions in each step. This argument is not used if \code{theta!=NULL} or \code{approach!=NULL}.
		%%     ~~Describe \code{H} here~~
	}
	\item{p}{
		The number of origin updates. This argument is not used if \code{theta!=NULL} or \code{approach!=NULL}.
		%%     ~~Describe \code{p} here~~
	}
	\item{thetaList}{
		A vector with the possible values for \code{theta}. This argument is not used if \code{theta} argument is not \code{NULL}.
		%%     ~~Describe \code{thetaList} here~~
	}
	\item{mc.cores}{
		Number of cores that will be used for estimate the theta parameter. It is not accepted \code{mc.cores>1} on Windows SO.
		%%     ~~Describe \code{mc.cores} here~~
	}
	\item{\dots}{
		Additional arguments for \code{tLineExtrap}.
		%%     ~~Describe \code{\dots} here~~
	}
}
\details{
	These functions are fully automatic, you just need to pass your time series. Particular cases are obtained by:
	If \code{theta = 1} the \code{tLineExtrapModel} method is computed;
	If \code{theta = 2} so the Standard Theta Method of Assimakopoulos and Nikolopoulos (2000) is computed.

	By default (\code{s=NULL}), the 90\% significance seasonal Z-test, used by Assimakopoulos and Nikolopoulos (2000), is applied for quarterly and monthly time series.
}
\value{
	An list containing the elements:
	\item{$y }{The original time series.}
    \item{$mean }{A time series element with the forecasting points.}
    \item{$fitted }{A time series element with the fitted points.}
    \item{$residuals }{A time series element with the residual points.}
    \item{$theta }{The estimated theta value.}
	\item{$tLineExtrap_par }{The estimated parameters of \code{tLineExtrap} method.}
    \item{$weights }{The estimated weights values.}
}
\references{
	Fioruci J.A., Pellegrini T.R., Louzada F., Petropoulos F. (2015). \emph{The Optimised Theta Method}. arXiv preprint, arXiv:1503.03529.

	Assimakopoulos, V. and Nikolopoulos k. (2000). \emph{The theta model: a decomposition approach to forecasting}. International Journal of Forecasting 16, 4, 521-530.
}
\author{
	Jose Augusto Fiorucci, Francisco Louzada
}
\note{
	The \code{thetaM} function is just a particular case of \code{otm} with \code{theta=2}.
}

	%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link[forecTheta]{forecTheta-package}}, \code{\link[forecTheta]{dotm}}, \code{\link[forecTheta]{groe}}
}
\examples{

y1 = 2+ 0.15*(1:20) + rnorm(20,2)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30,2)
y =  as.ts(c(y1,y2))

otm.arxiv(y, h=10)

### running the M3-competition data base by OTM approach (a) ###
#require(Mcomp)
#data(M3)
#
#forec = matrix(NA, nrow=3003, ncol=18)
#obs = matrix(NA, nrow=3003, ncol=18) #matrix of the out-sample values
#
#for(i in 1:3003){
#	if(i %% 100 == 0){print(i)}
#	x=M3[[i]]$x
#	h=M3[[i]]$h
#	out = otm.arxiv(x,h,approach='a',tLineExtrap=ses)
#	forec[i,1:h] = out$mean
#	obs[i,1:h] = M3[[i]]$xx
#}
#
#sAPE = errorMetric(obs, forec, type="sAPE", statistic="N") ## sAPE matrix
#
##### sMAPE results ##
### Yearly
#mean( sAPE[1:645, 1:6] )
### QUARTERLY
#mean( sAPE[646:1401, 1:8] )
### MONTHLY
#mean( sAPE[1402:2829, 1:18] )
### Other
#mean( sAPE[2830:3003, 1:8] )
### ALL
#mean( sAPE, na.rm=TRUE )

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{otm}
\keyword{thetaM}
\keyword{theta-method} % __ONLY ONE__ keyword per line
\keyword{time series forecasting} % __ONLY ONE__ keyword per line
