otm <- function( y, h=5, seasonal=NULL, theta=NULL, g="sAPE", n1=length(y)-h, 
						m=floor(h/2), H=h, p=1+floor((length(y)-n1)/m),
						thetaList=seq(from=1,to=5,by=0.5), tLineExtrap=ses, 
						mc.cores=1, ...){

	if(!is.ts(y)){ stop("ERROR in otm function: y must be a object of time series class.") }
	if(!is.numeric(h)){	stop("ERROR in otm function: h must be a positive integer number.")	}
	if(!is.null(theta) && !is.numeric(theta)){ stop("ERROR in otm function: theta must be a numeric value higher or equal to 1 or NULL.")}
	if(is.numeric(theta) && theta < 1){	stop("ERROR in otm function: theta must be a numeric value higher or equal to 1 or NULL.")}
	if(is.null(theta) && n1 < 4) stop("ERROR in otm function: n1 < 4)") 
	if(is.null(theta) && p > 1+floor((length(y)-n1)/m)){ stop("ERROR in otm function: p > 1+floor((length(y)-n1)/m)") }
	
	n = length(y)
	fq = frequency(y)
	time_y = time(y)
	
	if( is.null(seasonal) && any(fq == c(4,12)) ){
		xacf = acf(y,plot=FALSE)$acf[-1,1,1]
		clim = 1.64/sqrt(length(y))*sqrt(cumsum(c(1,2*xacf^2)))
		seasonal = abs(xacf[fq]) > clim[fq]
	}else{ 
		seasonal = FALSE 
	}
	
	if(seasonal){
        y_decomp = decompose(y, type = "multiplicative")$seasonal
        y = y/y_decomp
    }
	
	if(is.null(theta)){
		
		lossFunctionList = as.numeric(
			mclapply( X=thetaList, FUN=function(theta) groe( y=y, forecFunction=otm, g=g, 
														n1=n1, m=m, H=H, p=p, theta=theta, tLineExtrap=tLineExtrap, ...), 
						mc.cores = mc.cores )
		)
		
		aux = which.min(lossFunctionList)
		theta = thetaList[aux]
	}
	
	time_forec = time_y[n] + (1:h)/fq
	omega = 1 - 1/theta
	
	## linear method
	l = lm(y ~ time_y)
	l$mean = l$coeff[1] + l$coeff[2]*time_forec
	
	## other extrapolation method
	Z = l$fitted.values  +   theta * l$residuals
    X = tLineExtrap(Z, h=h, ...)

	Y_fitted = as.numeric( omega*l$fitted.values + (1-omega)*X$fitted )
	Y_fitted = ts(Y_fitted, start = start(y), frequency = frequency(y))
    Y_fcast = as.numeric( omega*l$mean + (1-omega)*X$mean )
	Y_fcast = ts(Y_fcast, start = end(y)+c(0,1), frequency = frequency(y))
    Y_residuals = y - Y_fitted
	
	if(seasonal){
		y = y*y_decomp
        Y_fitted = y_decomp * Y_fitted
        Y_fcast = snaive(y_decomp, h = h)$mean * Y_fcast
        Y_residuals = y - Y_fitted
    }
	
    fit = list()
	fit$seasTest = seasonal
    fit$mean = Y_fcast
    fit$fitted = Y_fitted
    fit$residuals = Y_residuals
	fit$theta = theta
	fit$weights = c(omega, 1-omega)
	
	return(fit)
}




thetaM <- function( y, h=5, seasonal=NULL){
	if(!is.ts(y)){ stop("ERROR in thetaM function: y must be a object of time series class.") }
	if(!is.numeric(h)){stop("ERROR in thetaM function: h must be a positive integer number.")}
	
	fit = otm( y=y, h=h, seasonal=seasonal, theta=2)
	
	return( fit )
}

