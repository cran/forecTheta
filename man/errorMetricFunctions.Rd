\name{errorMetric}
\alias{errorMetric}
	%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Error Metric Function
	%%  ~~function to do ... ~~
}
\description{
	This function implements some of the more used error metrics. 
	These metrics are "sMAPE", "MAPE", "MAE", "MSE" and they respectively versions with median  
	"sMdAPE", "MdAPE", "MdAE", "MdSE".
	%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
	errorMetric(obs, forec, type="sAPE", statistic="M")
}
	%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
	A vector or a matrix with the real values.
	%%     ~~Describe \code{obs} here~~
	}
  \item{forec}{
	A vector or a matrix with the estimated values.
	%%     ~~Describe \code{forec} here~~
}
  \item{type}{
	The error type of "sAPE", "APE", "AE" and "SE".
	%%     ~~Describe \code{type} here~~
}
  \item{statistic}{
	The statistic to be returned. Use "M" or "Md" for return the mean or median of the errors. 
	If "N" so a vector with all errors will be returned.  
	%%     ~~Describe \code{statistic} here~~
}
}
\details{
	The metric sMAPE is obtained using \code{type = "sAPE"} and  \code{statistic = "M"}
	
	The metric sMdAPE is obtained using \code{type = "sAPE"} and  \code{statistic = "Md"}
	
	The metric MAPE is obtained using \code{type = "APE"} and  \code{statistic = "M"}
	
	The metric MdAPE is obtained using \code{type = "APE"} and  \code{statistic = "Md"}
	
	The metric MAE is obtained using \code{type = "AE"} and  \code{statistic = "M"}
	
	The metric MdAE is obtained using code{type = "AE"} and  \code{statistic = "Md"}
	
	The metric MSE is obtained using \code{type = "SE"} and  \code{statistic = "M"}
	
	The metric MdSE is obtained using \code{type = "SE"} and  \code{statistic = "Md"}
	%%  ~~ If necessary, more details than the description above ~~
}
\value{
	If \code{statistic="M"} or \code{statistic="Md"} it is returned the respectively error metric result.
	If \code{statistic="N"} so is returned a vector with all errors points according to the chosen error type.  
	%%  ~Describe the value returned
	%%  If it is a LIST, use
	%%  \item{comp1 }{Description of 'comp1'}
	%%  \item{comp2 }{Description of 'comp2'}
	%% ...
}
%%\references{
	%% ~put references to the literature/web site here ~
%%}
\author{
	Jose Augusto Fioruci and Francisco Louzada
	%%  ~~who you are~~
}
%%\note{
	%%  ~~further notes~~
%%}

	%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
	\code{\link[forecTheta]{forecTheta-package}}
	%%~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##############################################################	

y1 = 2+ 0.15*(1:20) + rnorm(20,2)
y2 = y1[20]+ 0.3*(1:30) + rnorm(30,2)
y =  as.ts(c(y1,y2))

otm.fit <- otm(y=as.ts(y[1:40]), h=10)

### sMAPE metric
errorMetric(obs=as.ts(y[41:50]), forec=otm.fit$mean)

### sMdAPE metric
errorMetric(obs=as.ts(y[41:50]), forec=otm.fit$mean, statistic = "Md")

### MASE metric
meanDiff1 = mean(abs(diff(as.ts(y[1:40]), lag = 1)))
errorMetric(obs=as.ts(y[41:50]), forec=otm.fit$mean, type = "AE", statistic = "M") / meanDiff1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ error metric }
\keyword{ sMAPE }
\keyword{ sMdAPE }
\keyword{ MAE }
\keyword{ MdAE }
\keyword{ MSE }
\keyword{ MdSE }% __ONLY ONE__ keyword per line
