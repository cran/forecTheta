\name{PI_eval}
\alias{PI_eval}
\title{Prediction Interval (PI) evaluation}
\description{
Implementation of Mean Scaled Interval Score (MSIS) and Absolute Coverage Difference (ACD).
}
\usage{
PI_eval(obs, forec, lower_bounds, upper_bounds, name = c("MSIS", "ACD"), alpha = 0.05)
}
\arguments{
  \item{obs}{
    Observed time series values.
  }
  \item{forec}{
    Forecasted values corresponding to the observed data.
  }
  \item{lower_bounds}{
    Lower bounds of the prediction interval.
  }
  \item{upper_bounds}{
    Upper bounds of the prediction interval.
  }
  \item{name}{
    Evaluation metric to be used. Options are \code{"MSIS"} or \code{"ACD"}.
  }
  \item{alpha}{
    Significance level for the prediction interval. Default is 0.05.
  }
}
\details{
The MSIS is computed by combining the width of the prediction interval with a penalty for
observations falling outside the interval. This score is averaged over time and scaled by
the mean absolute seasonal difference, following the same approach used for the MASE in
the M4 competition, to ensure scale invariance. As a complementary metric, the ACD
measures the absolute difference between the empirical and nominal coverage (e.g., 0.95),
helping to assess interval reliability.
}
\value{
Returns a numeric value corresponding to the selected evaluation metric (MSIS or ACD).
}
\references{
Gneiting, T. and Raftery, A.E. (2007). Strictly proper scoring rules, prediction, and estimation. Journal of the American Statistical Association 102 (477), 359--378, \doi{10.1198/016214506000001437}.

Makridakis, S., Spiliotis, E. and Assimakopoulos, V. (2018). The M4 Competition: 100,000 time series and 61 forecasting methods. International Journal of Forecasting 36 (1), 54--74, \doi{10.1016/j.ijforecast.2019.04.014}.
}
\author{
Jose Augusto Fiorucci, Francisco Louzada, Igor De Oliveira Barros Faluhelyi.
}
\seealso{
\code{\link[forecTheta]{forecTheta-package}},
\code{\link[forecTheta]{errorMetricFunctions}},
\code{\link[forecTheta]{dotm}},
\code{\link[forecTheta]{dstm}},
\code{\link[forecTheta]{otm}},
\code{\link[forecTheta]{stm}}
}
\examples{
y1 <- 2 + 0.15 * (1:20) + rnorm(20)
y2 <- y1[20] + 0.3 * (1:30) + rnorm(30)
y <- as.ts(c(y1, y2))
out <- dotm(y, h = 10)

# MSIS metric
PI_eval(
  obs = y, forec = out$mean,
  lower_bounds = out$lower[, 3], upper_bounds = out$upper[, 3],
  name = "MSIS", alpha = 0.05
)

# ACD metric
PI_eval(
  obs = y, forec = out$mean,
  lower_bounds = out$lower[, 3], upper_bounds = out$upper[, 3],
  name = "ACD", alpha = 0.05
)
}
\keyword{MSIS}
\keyword{ACD}
\keyword{time series forecasting}
