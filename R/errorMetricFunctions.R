errorMetric <- function(obs, forec, type="sAPE", statistic="M"){

	if( !any(type==c("AE","SE","APE","sAPE")) ) stop("Error in errorMetric function: this error type has not been implemented.")
	if( !any(statistic==c("M","Md","N")) ) stop("Error in errorMetric function: this statistic has not been implemented.")
	if( is.ts(obs) || is.ts(forec) ){
		obs = as.numeric(obs)
		forec = as.numeric(forec)
	}
	obs = as.matrix(obs)
	forec = as.matrix(forec)
	if(any(dim(obs)!=dim(forec)))  stop("Error in errorMetric function: the dimensions of the vectors are different.")

	if(type == "AE")
		errors = abs(obs - forec)

	if(type == "SE")
		errors = (obs - forec)^2

	if(type == "APE")
		errors = abs( 100*(obs - forec)/obs )

	if(type == "sAPE")
		errors =  abs( 200*(obs - forec)/(abs(obs) + abs(forec)) )

	if(statistic == "M")
		return(  mean(errors, na.rm=TRUE) )

	if(statistic == "Md")
		return(  median(errors, na.rm=TRUE) )

	return( errors )

}

#########################################################################################################################

PI_eval <- function(obs, forec, lower_bounds, upper_bounds, name = c("MSIS", "ACD"), alpha = .05){

  if( !any(name==c("MSIS", "ACD")) ) stop("Error in Prediction Interval (PI) evaluation function: this error metric is not implemented.")
  if( !any( (alpha < 1) && (alpha > 0) ) ) stop("Error in Prediction Interval (PI) evaluation function: alpha parameter should be greater than 0 and less than 1")

  obs <- if (is.ts(obs)) obs else as.ts(obs)
  forec <- if (is.ts(forec)) forec else as.ts(forec)
  lower_bounds <- if (is.ts(lower_bounds)) lower_bounds else as.ts(lower_bounds)
  upper_bounds <- if (is.ts(upper_bounds)) upper_bounds else as.ts(upper_bounds)

  if(name == "MSIS"){

    snaive_forecast <- snaive(obs, h = length(forec))$mean
    scale <- mean(
      abs(
        as.numeric(snaive_forecast) - as.numeric(forec)
      )
    )

    interval_score <- mean(
      as.numeric((upper_bounds - lower_bounds)) +
        (2 / alpha) * as.numeric((lower_bounds - forec)) * as.numeric((forec < lower_bounds)) +
        (2 / alpha) * as.numeric((forec - upper_bounds)) * as.numeric((forec > upper_bounds))
    )

    value = interval_score / scale

  }

  if(name == "ACD"){
    coverage_rate = mean(
      as.numeric((forec >= lower_bounds)) & as.numeric((forec <= upper_bounds))
    )

    value = abs(
      coverage_rate - (1 - alpha)
    )

  }

  return( value )

}
