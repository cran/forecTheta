\name{Bagged Theta Models}
\alias{bagged_dotm}
\alias{bagged_dstm}
\alias{bagged_otm}
\alias{bagged_stm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Bagged Theta Models
	%%  ~~function to do ... ~~
}
\description{
	%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Bagged implementation (Bergmeir et al, 2016) of Dynamic Optimised Theta Model,
Dynamic Standard Theta Model, Optimised Theta Model and
Standard Theta Model (Fiorucci et al, 2016).
}
\usage{
bagged_dotm(y, h=5, level=c(80,90,95), num_bootstrap = 100, bs_bootstrap = NULL,
  s_type="multiplicative", s_test="default", lambda=NULL,
  par_ini=c(y[1]/2, 0.5, 2),	estimation=TRUE, lower=c(-1e+10, 0.1, 1.0),
  upper=c(1e+10, 0.99, 1e+10), opt.method="Nelder-Mead", xreg=NULL)

bagged_dstm(y, h=5, level=c(80,90,95), num_bootstrap = 100, bs_bootstrap = NULL,
  s_type="multiplicative", s_test="default", lambda=NULL, par_ini=c(y[1]/2, 0.5),
  estimation=TRUE, lower=c(-1e+10, 0.1), upper=c(1e+10, 0.99),
  opt.method="Nelder-Mead", xreg=NULL)

bagged_otm(y, h=5, level=c(80,90,95), num_bootstrap = 100, bs_bootstrap = NULL,
  s_type="multiplicative", s_test="default", lambda=NULL, par_ini=c(y[1]/2, 0.5, 2),
  estimation=TRUE, lower=c(-1e+10, 0.1, 1.0), upper=c(1e+10, 0.99, 1e+10),
  opt.method="Nelder-Mead", xreg=NULL)

bagged_stm(y, h=5, level=c(80,90,95), num_bootstrap = 100, bs_bootstrap = NULL,
  s_type="multiplicative", s_test="default", lambda=NULL, par_ini=c(y[1]/2, 0.5),
  estimation=TRUE, lower=c(-1e+10, 0.1), upper=c(1e+10, 0.99),
  opt.method="Nelder-Mead", xreg=NULL)
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
	\item{num_bootstrap}{
    Number of bootstrap samples to be considered in the forecast (\code{num} argument of the \code{forecast::bld.mbb.bootstrap} function).
	}
	\item{bs_bootstrap}{
    Block size for bootstrap samples to be considered in the forecast (\code{block_size} argument of the \code{forecast::bld.mbb.bootstrap} function). Used only if \code{num_bootstrap >= 1}.
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

}
\details{
This version allows bagging to be used in conjunction with the models DOTM, DSTM, OTM and STM (Fiorucci et al, 2016), which usually generates considerably more accurate results, whether point or interval forecasts. In this case, \emph{Box-Cox and Loess-based decomposition bootstrap (BLD-MBB)} is used, for details see Bergmeir et al, 2016.

If \code{num_bootstrap > 1}, then the bagged series are extrapolated through simulation. In the end, point forecasts are calculated as the mean (\code{$mean}) or the median (\code{$median}) of these series, while interval forecasts are computed based on quantiles. To ensure a substantial number of simulations, each bagged series can be extrapolated multiple times through simulation, aiming for a total simulation volume of at least 10,000 series.

By default, the 90\% significance seasonal Z-test, used by Assimakopoulos and Nikolopoulos (2000), is applied for quarterly and monthly time series. The possibility of first checking the unit root was included because it was pointed out that this test presents many flaws for time series with this characteristic (Fiorucci et al, 2016). In this case, the KPSS test is performed with a significance level of 5\% and in the case of a unit root, then the series is differentiated before checking for seasonal behavior.

}
\value{
	An object of \code{thetaModel} class with one list containing the elements:
	\item{$method}{The name of the model/method}
	\item{$y}{The original time series.}
	\item{$s}{A binary indication for seasonal decomposition (original time series).}
	\item{type}{Seasonal decomposition type (original time series).}
	\item{opt.method}{The optimisation method used in the \code{optim()} function.}
  \item{$par }{The estimated values of \code{(ell, alpha, theta)} parameters (original time series)}
	\item{$weights }{The estimated weights values (original time series).}
  \item{$fitted }{A time series element with the fitted points (original time series).}
  \item{$residuals }{A time series element with the residual points (original time series).}
  \item{$mean}{The forecasting values calculed as the mean of bootstrapping simulations.}
  \item{$median}{The forecasting values calculed as the median of bootstrapping simulations.}
	\item{$level}{The levels for prediction intervals.}
	\item{$lower}{Lower limits for prediction intervals.}
	\item{$upper}{Upper limits for prediction intervals.}
}

\references{
Assimakopoulos, V. and Nikolopoulos k. (2000). \emph{The theta model: a decomposition approach to forecasting}. International Journal of Forecasting 16 (4), 521--530, <doi:10.1016/S0169-2070(00)00066-2>.

Bergmeir, C., Hyndman, R.J. and Benítez, J. M. (2016). \emph{Bagging exponential smoothing methods using STL decomposition and Box–Cox transformation}. International journal of forecasting 32 (2), 303--312, <doi:10.1016/j.ijforecast.2015.07.002>.

Fiorucci J.A., Pellegrini T.R., Louzada F., Petropoulos F., Koehler, A. (2016). \emph{Models for optimising the theta method and their relationship to state space models},	International Journal of Forecasting, 32 (4), 1151--1161, <doi:10.1016/j.ijforecast.2016.02.005>.
}

\author{
	Jose Augusto Fiorucci, Francisco Louzada, Igor De Oliveira Barros Faluhelyi.
}

\seealso{
\code{\link[forecTheta]{forecTheta-package}}, \code{\link[forecTheta]{dotm}},
\code{\link[forecTheta]{dstm}}, \code{\link[forecTheta]{otm}},
\code{\link[forecTheta]{stm}}
}

\examples{


#### additive seasonal decomposition ###
x = sin(2*pi*seq(0,9,len=300)) + exp((1:300)/150) + rnorm(mean=100,sd=0.5,n=300)
y = ts(x, frequency=33)
out <- bagged_dotm(y, h=50, s_type='additive', num_bootstrap = 5)
summary(out)
plot(out)


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DOTM}
\keyword{DSTM}
\keyword{OTM}
\keyword{STM}
\keyword{time series forecasting} % __ONLY ONE__ keyword per line
