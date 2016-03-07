\name{Theta Models}
\alias{dotm}
\alias{dstm}
\alias{otm}
\alias{stm}
\alias{stheta}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Theta Models
	%%  ~~function to do ... ~~
}
\description{
	%%  ~~ A concise (1-5 lines) description of what the function does. ~~
	Functions for forecast univariate time series using the Dynamic Optimised Theta Model, Dynamic Standard Theta Model,
	Optimised Theta Model and Standard Theta Model (Fiorucci et al, 2016). 
	We also provide an implementation for the Standard Theta Method (STheta) of Assimakopoulos and Nikolopoulos (2000).
}
\usage{
	dotm(y, h=5, level=c(80,90,95), s=NULL, par_ini=c(y[1]/2, 0.5, 2), 
	estimation=TRUE, lower=c(-1e+10, 0.1, 1.0), upper=c(1e+10, 0.99, 1e+10))
	
	dstm(y, h=5, level=c(80,90,95), s=NULL, par_ini=c(y[1]/2, 0.5), 
	estimation=TRUE, lower=c(-1e+10, 0.1), upper=c(1e+10, 0.99))
	
	otm(y, h=5, level=c(80,90,95), s=NULL, par_ini=c(y[1]/2, 0.5, 2), 
	estimation=TRUE, lower=c(-1e+10, 0.1, 1.0), upper=c(1e+10, 0.99, 1e+10))
	
	stm(y, h=5, level=c(80,90,95), s=NULL, par_ini=c(y[1]/2, 0.5), 
	estimation=TRUE, lower=c(-1e+10, 0.1), upper=c(1e+10, 0.99))
	
	stheta(y, h=5, s=NULL)
}

\arguments{
	\item{y}{
		Object of time series class.
	}
	\item{h}{
		Number of required forecasting periods.
	}
	\item{level}{
		Levels for prediction intervals.
	}
	\item{s}{
		If \code{TRUE}, the multiplicative seasonal decomposition is used. 
		If \code{NULL}, quarterly and monthly time series are tested for statistically seasonal behaviour, with 90\% of significance. Default is NULL.
	}
	\item{par_ini}{
		Vector of initialization for \code{(ell, alpha, theta)} parameters.
	}
	\item{estimation}{
		If \code{TRUE}, the \code{optim()} function is consider for compute the minimum square estimator of parameters.
		If \code{FALSE}, the models/methods are computed for \code{par_ini} values.
	}
	\item{lower}{
		The lower limit of parametric space.
	}
	\item{upper}{
		The upper limit of parametric space.
	}
	
}
\details{
	By default (\code{s=NULL}), the 90\% significance seasonal Z-test, used by Assimakopoulos and Nikolopoulos (2000), is applied for quarterly and monthly time series. 

	For details of each model see Fiorucci et al, 2016. 
	If you are looking for the methods presented in the arXiv paper (Fiorucci et al, 2015), see \code{otm.arxiv()} function.
	}
\value{
	An object of \code{thetaModel} class with one list containing the elements:
	\item{$method}{The name of the model/method}
	\item{$y }{The original time series.}
	\item{$s}{A binary indication for seasonal decomposition}
    \item{$par }{The estimated values for \code{(ell, alpha, theta)} parameters}
	\item{$weights }{The estimated weights values.}
    \item{$fitted }{A time series element with the fitted points.}
    \item{$residuals }{A time series element with the residual points.}
    \item{$mean}{The forecasting values}
	\item{$level}{The levels for prediction intervals}
	\item{$lower}{Lower limits for prediction intervals}
	\item{$upper}{Upper limits for prediction intervals}  
}
\references{
	Fiorucci J.A., Pellegrini T.R., Louzada F., Petropoulos F., Koehler, A. (2016). \emph{Models for optimising the theta method and their relationship to state space models}, 
	International Journal of Forecasting. Accepted Paper. \url{https://www.researchgate.net/publication/294420765_Models_for_optimising_the_theta_method_and_their_relationship_to_state_space_models}
	
	Fiorucci J.A., Pellegrini T.R., Louzada F., Petropoulos F. (2015). \emph{The Optimised Theta Method}. Free available at \url{ http://arxiv.org/abs/1503.03529 }. 
	
	Assimakopoulos, V. and Nikolopoulos k. (2000). \emph{The theta model: a decomposition approach to forecasting}. International Journal of Forecasting 16, 4, 521-530.	
}
\author{
	Jose Augusto Fiorucci, Francisco Louzada and Bao Yiqi
}

\seealso{
	\code{\link[forecTheta]{forecTheta-package}}, \code{\link[forecTheta]{otm.arxiv}}
}

\examples{

y1 = 2+ 0.15*(1:20) + rnorm(20,2)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30,2)
y =  as.ts(c(y1,y2))

out <- dotm(y, h=10)
summary(out)
plot(out)

######### Reproducing the M3 results by DOTM ############
# library(Mcomp)
# data(M3)
#
# forec = matrix(NA, nrow=3003, ncol=18)
# obs = matrix(NA, nrow=3003, ncol=18) #matrix of the out-sample values
# meanDiff <- rep(1, 3003)
# 
# for(i in 1:3003){
#	 if(i \%\% 100 == 0){print(i);}
#	 x=M3[[i]]$x
#	 h=M3[[i]]$h
#	 out = dotm(x,h)
#	 forec[i,1:h] = out$mean
#	 obs[i,1:h] = M3[[i]]$xx
#	 meanDiff[i] = mean(abs(diff(x, lag = frequency(x))))
# }

############## sMAPE ###################
#	sAPE_matrix = errorMetric(obs=obs, forec=forec, type="sAPE", statistic="N")
#### Yearly ###
#	mean( sAPE_matrix[1:645, 1:6] )
#### QUARTERLY ###
#	mean( sAPE_matrix[646:1401, 1:8] )
#### MONTHLY ###
#	mean( sAPE_matrix[1402:2829, 1:18] )
#### Other ###
#	mean( sAPE_matrix[2830:3003, 1:8] )
#### ALL ###
#	mean( sAPE_matrix, na.rm=TRUE )
#	
############# MASE ######################	
#	AE_matrix = errorMetric(obs=obs, forec=forec, type="AE", statistic="N")
#   ASE_matrix=AE_matrix/meanDiff
#### Yearly ###
#	mean( ASE_matrix[1:645, 1:6] )
#### QUARTERLY ###
#	mean( ASE_matrix[646:1401, 1:8] )
#### MONTHLY ###
#	mean( ASE_matrix[1402:2829, 1:18] )
#### Other ###
#	mean( ASE_matrix[2830:3003, 1:8] )
#### ALL ###
#	mean( ASE_matrix, na.rm=TRUE )
########################################################	
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DOTM}
\keyword{DSTM}
\keyword{OTM}
\keyword{STM}
\keyword{STheta} % __ONLY ONE__ keyword per line
\keyword{time series forecasting} % __ONLY ONE__ keyword per line
