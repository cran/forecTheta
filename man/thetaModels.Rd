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
dotm(y, h=5, level=c(80,90,95),	s_type="multiplicative", s_test="default",
  	lambda=NULL, par_ini=c(y[1]/2, 0.5, 2),	estimation=TRUE,
  	lower=c(-1e+10, 0.1, 1.0), upper=c(1e+10, 0.99, 1e+10),
  	opt.method="Nelder-Mead", xreg=NULL, s=NULL)

dstm(y, h=5, level=c(80,90,95), s_type="multiplicative", s_test="default",
  	lambda=NULL, par_ini=c(y[1]/2, 0.5), estimation=TRUE,
  	lower=c(-1e+10, 0.1), upper=c(1e+10, 0.99),
  	opt.method="Nelder-Mead", xreg=NULL, s=NULL)

otm(y, h=5, level=c(80,90,95), s_type="multiplicative", s_test="default",
  	lambda=NULL, par_ini=c(y[1]/2, 0.5, 2), estimation=TRUE,
  	lower=c(-1e+10, 0.1, 1.0), upper=c(1e+10, 0.99, 1e+10),
  	opt.method="Nelder-Mead", xreg=NULL, s=NULL)

stm(y, h=5, level=c(80,90,95), s_type="multiplicative", s_test="default",
  	lambda=NULL, par_ini=c(y[1]/2, 0.5), estimation=TRUE,
  	lower=c(-1e+10, 0.1), upper=c(1e+10, 0.99),
  	opt.method="Nelder-Mead", xreg=NULL, s=NULL)

stheta(y, h=5, s_type="multiplicative", s_test="default", s=NULL)
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
	\item{s_type}{
	  Use \code{"multiplicative"} for the classic multiplicative seasonal decomposition,
	  \code{"additive"} for the classic additive seasonal decomposition and
	  \code{"stl"} for the STL decomposition (in this case, \code{forecast::mstl} is considered).
	}
	\item{s_test}{
		If \code{TRUE}, seasonal decomposition is used.
		If \code{FALSE}, seasonal decomposition is not used.
		If \code{default}, the time series is tested for statistically seasonal behaviour.
		If \code{uni_root}, then first a difference is taken if a unit root is detected.
	}
	\item{lambda}{
	Parameter for Box-Cox transformation (\code{forecast::BoxCox} is considered).
	If \code{lambda=NULL}, then it is ignored.
	If \code{lambda="auto"}, then this parameter is automatically selected (see \code{forecast::BoxCox.lambda()}).
	Furthermore, \code{lambda=0} corresponds to the logarithmic transformation, which can be used to get strictly positive forecasts.
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
	\item{opt.method}{
		The numeric optimisation method for \code{optim()} function.
		Choose one among \code{'Nelder-Mead', 'L-BFGS-B', 'SANN'}.
	}
	\item{xreg}{
		A matrix with the regressor variables including the out-of-sample data.
	}
	\item{s}{
		The argument \code{s} is deprecated and has been replaced by \code{s_type} and \code{s_test} for improved functionality and clarity. While \code{s} is still supported in this version for backward compatibility, it may be removed in future releases. Users are strongly encouraged to update their code accordingly.
	}

}
\details{
  By default, the 90\% significance seasonal Z-test, used by Assimakopoulos and Nikolopoulos (2000), is applied for quarterly and monthly time series. The possibility of first checking the unit root was included because it was pointed out that this test presents many flaws for time series with this characteristic (Fiorucci et al, 2016). In this case, the KPSS test is performed with a significance level of 5\% and in the case of a unit root, then the series is differentiated before checking for seasonal behavior.

For details of each model see Fiorucci et al, 2016.
If you are looking for the methods presented in the arXiv paper (Fiorucci et al, 2015), see \code{\link[forecTheta]{otm.arxiv}} function.

This version allows bagging to be used in conjunction with these models, see \code{\link[forecTheta]{bagged_dotm}},
	\code{\link[forecTheta]{bagged_dstm}},
	\code{\link[forecTheta]{bagged_otm}} and
	\code{\link[forecTheta]{bagged_stm}} functions.
}
\value{
	An object of \code{thetaModel} class with one list containing the elements:
	\item{$method}{The name of the model/method}
	\item{$y }{The original time series.}
	\item{$s}{A binary indication for seasonal decomposition.}
	\item{type}{Classical seasonal decomposition type.}
	\item{opt.method}{The optimisation method used in the \code{optim()} function.}
    \item{$par }{The estimated values for \code{(ell, alpha, theta)} parameters}
	\item{$weights }{The estimated weights values.}
    \item{$fitted }{A time series element with the fitted points.}
    \item{$residuals }{A time series element with the residual points.}
    \item{$mean}{The forecasting values.}
	\item{$level}{The levels for prediction intervals.}
	\item{$lower}{Lower limits for prediction intervals.}
	\item{$upper}{Upper limits for prediction intervals.}
	\item{$tests}{The p.value of Teraesvirta Neural Network test applied on unseasoned time series and the p.value of Shapiro-Wilk test applied on unseasoned residuals.}
}
\references{
	Fiorucci J.A., Pellegrini T.R., Louzada F., Petropoulos F., Koehler, A. (2016). \emph{Models for optimising the theta method and their relationship to state space models},	International Journal of Forecasting, 32 (4), 1151--1161, <doi:10.1016/j.ijforecast.2016.02.005>.

	Assimakopoulos, V. and Nikolopoulos k. (2000). \emph{The theta model: a decomposition approach to forecasting}. International Journal of Forecasting 16 (4), 521--530, <doi:10.1016/S0169-2070(00)00066-2>.

}
\author{
	Jose Augusto Fiorucci, Francisco Louzada
}

\seealso{
	\code{\link[forecTheta]{forecTheta-package}},
	\code{\link[forecTheta]{bagged_dotm}},
	\code{\link[forecTheta]{bagged_dstm}},
	\code{\link[forecTheta]{bagged_otm}},
	\code{\link[forecTheta]{bagged_stm}},
	\code{\link[forecTheta]{otm.arxiv}}
}

\examples{

y1 = 2+ 0.15*(1:20) + rnorm(20)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30)
y =  as.ts(c(y1,y2))
out <- dotm(y, h=10)
summary(out)
plot(out)

#### additive seasonal decomposition ###
x = sin(2*pi*seq(0,9,len=300)) + exp((1:300)/150) + rnorm(mean=0,sd=0.5,n=300)
y = ts(x, frequency=33)
out <- dotm(y, h=50, s_type='additive')
summary(out)
plot(out)

# #########################################################
# ######### Reproducing the M3 results by DOTM ############
# #########################################################
#
# library(Mcomp)
# data(M3)
#
# forec = matrix(NA, nrow=3003, ncol=18)
# obs = matrix(NA, nrow=3003, ncol=18) #matrix of the out-sample values
# meanDiff <- rep(1, 3003)
#
# for(i in 1:3003){
#   x=M3[[i]]$x
#   h=M3[[i]]$h
#   out = dotm(x,h,level=NULL)
#   forec[i,1:h] = out$mean
#   obs[i,1:h] = M3[[i]]$xx
#   meanDiff[i] = mean(abs(diff(x, lag = frequency(x))))
# }
#
# ############## sMAPE ###################
# sAPE_matrix = errorMetric(obs=obs, forec=forec, type="sAPE", statistic="N")
# #### Yearly ###
# mean( sAPE_matrix[1:645, 1:6] )
# #### QUARTERLY ###
# mean( sAPE_matrix[646:1401, 1:8] )
# #### MONTHLY ###
# mean( sAPE_matrix[1402:2829, 1:18] )
# #### Other ###
# mean( sAPE_matrix[2830:3003, 1:8] )
# #### ALL ###
# mean( sAPE_matrix, na.rm=TRUE )
# #
# ############# MASE ######################
# AE_matrix = errorMetric(obs=obs, forec=forec, type="AE", statistic="N")
# ASE_matrix=AE_matrix/meanDiff
# #### Yearly ###
# mean( ASE_matrix[1:645, 1:6] )
# #### QUARTERLY ###
# mean( ASE_matrix[646:1401, 1:8] )
# #### MONTHLY ###
# mean( ASE_matrix[1402:2829, 1:18] )
# #### Other ###
# mean( ASE_matrix[2830:3003, 1:8] )
# #### ALL ###
# mean( ASE_matrix, na.rm=TRUE )
# ########################################################


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DOTM}
\keyword{DSTM}
\keyword{OTM}
\keyword{STM}
\keyword{STheta} % __ONLY ONE__ keyword per line
\keyword{time series forecasting} % __ONLY ONE__ keyword per line
