\name{crossValidation}
\alias{groe}
\alias{fixOrig}
\alias{rolOrig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Generalised Rolling Origin Evaluation
}
\description{
	This function implements the Generalised Rolling Origin Evaluation of Fioruci et al (2015). 
	Its particular cases include the cross validation methods: Rolling Origin Evaluation and Fixed Origin Evaluation of Tashman(2000).
	%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
	groe(y, forecFunction=thetaM, g="sAPE", n1=length(y)-10, m=5, 
	        H=length(y)-n1, p=1+floor((length(y)-n1)/m), ...)
	
	rolOrig(y, forecFunction=thetaM, g="sAPE", n1=length(y)-10, ...)
	
	fixOrig(y, forecFunction=thetaM, g="sAPE", n1=length(y)-10, ...)
}
	%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{y}{
	Object of time series class or a vector
	%%     ~~Describe \code{y} here~~
	}
	\item{forecFunction}{
	A forecasting method as one object of the \code{forecast} class of forecast package. 
	Default is the function \code{thetaM}, which implements the Theta-method.
	%%     ~~Describe \code{forecFunction} here~~
	}
	\item{g}{
	The prediction error type of \code{errorMetric} function. The possible values are \code{"sAPE", "APE", "AE"} and \code{"SE"}.
	%%     ~~Describe \code{g} here~~
	}
	\item{n1}{
	The index of the first origin element. 
	%%     ~~Describe \code{n1} here~~
	}
	\item{m}{
	The number of movements of the origin in each update.
	%%     ~~Describe \code{m} here~~
	}
	\item{H}{
	The number of predictions forward of each origin.
	%%     ~~Describe \code{H} here~~
	}
	\item{p}{
		The number of origin updates. 
		Default is the maximum.
	%%     ~~Describe \code{p} here~~
	}
	\item{\dots}{
	Additional arguments for \code{forecFunction}.
	%%     ~~Describe \code{\dots} here~~
	}
}
\details{
	If \code{m=1} is computed the Rolling Origin Evaluation.
	If \code{m>=length(y)-n1} is computed the Fixed Origin Evaluation, where \code{h} is the required number of forecasting.
	%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The sum of the prediction errors.
	%%  ~Describe the value returned
	%%  If it is a LIST, use
	%%  \item{comp1 }{Description of 'comp1'}
	%%  \item{comp2 }{Description of 'comp2'}
	%% ...
}
\references{
	Fioruci J.A., Pellegrini T.R., Louzada F., Petropoulos F. (2015). \emph{The Optimised Theta Method}. Free available at \url{ http://arxiv.org/abs/1503.03529 }. 
	
	Tashman, L.J. (2000). \emph{Out-of-sample tests of forecasting accuracy: an analysis and review}. International Journal of Forecasting 16 (4), 437--450.
%% ~put references to the literature/web site here ~
}
\author{
	Jose Augusto Fioruci and Francisco Louzada
	%%  ~~who you are~~
}
\note{
	The \code{otm} function use this function for estimate the theta parameter when the \code{theta} argument is \code{NULL}. 
	Your computer may go into an infinite looping if you use \code{forecFunction = otm} and not specific a numeric value for the \code{theta} argument.
	%%  ~~further notes~~
}

	%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link[forecTheta]{forecTheta-package}}, \code{\link[forecTheta]{otm}}, \code{\link[forecast]{ses}} 
	%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
y1 = 2+ 0.15*(1:20) + rnorm(20,2)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30,2)
y =  as.ts(c(y1,y2))

## Rolling Origin Evaluation
rolOrig( y=y, forecFunction = ses, n1=40)
rolOrig( y=y, forecFunction = thetaM, n1=40)
rolOrig( y=y, forecFunction = otm, n1=40, theta=1)
rolOrig( y=y, forecFunction = otm, n1=40, theta=2)
rolOrig( y=y, forecFunction = otm, n1=40, theta=3)

## Fixed Origin Evaluation
fixOrig( y=y, forecFunction = ses, n1=40)
fixOrig( y=y, forecFunction = thetaM, n1=40)
fixOrig( y=y, forecFunction = otm, n1=40, theta=1)
fixOrig( y=y, forecFunction = otm, n1=40, theta=2)
fixOrig( y=y, forecFunction = otm, n1=40, theta=3)

## Generalised Rolling Origin Evaluation with two origin updates.
## Where the first is the 40th element and second is the 45th element
groe( y=y, forecFunction = ses, m=5, n1=40)
groe( y=y, forecFunction = thetaM, m=5, n1=40)
groe( y=y, forecFunction = otm, m=5, n1=40, theta=1)
groe( y=y, forecFunction = otm, m=5, n1=40, theta=2)
groe( y=y, forecFunction = otm, m=5, n1=40, theta=3)
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Generalised Rolling Origin Evaluation }
\keyword{ Fixed Origin Evaluation }
\keyword{ Rolling Origin Evaluation }
\keyword{ Cross Validation }% __ONLY ONE__ keyword per line