rm(list=ls())

library(PointFore)
library(isodistrreg)
library(tidyverse)
library(reshape2)

setClass("IDR_quantiles", slots=list(X="numeric", Y="numeric", quantiles="matrix", alpha="numeric"))
setOldClass('pointfore')
setClass("estimated_levels", slots=list(x.grid="numeric", alpha_x="numeric", alpha_gmm="numeric", 
                                        gmm="logical", res_parametric='pointfore', forec="numeric"))

#X: Vector of forecasts, Y: Vector of observations, 
#alpha: vector of quantile levels
estimate_quantiles <- function(forec, obs, alpha=seq(0.1,0.9,0.1)) {
  #Initialize vector for alpha(x)
  alpha_x <- numeric()
   
  #Fit IDR model
  fit <- idr(
    y = obs,
    X = data.frame(X=forec)
  )
  prediction <- predict(fit, data = data.frame(X=forec))
  quantiles <- qpred(prediction, quantiles = alpha)
  
  res <- new('IDR_quantiles',X=forec, Y=obs, quantiles=quantiles, alpha=alpha)
  
  return(res)
}

#x object of class IDR_quantiles (result of estimate_quantiles)
plot.quantiles <- function(x) {
  data <- data.frame('X'=x@X, 'Y'=x@Y)
  data <- cbind(data, x@quantiles)
  data <- melt(data, id.vars = c('X', 'Y'))
  p <- ggplot() + 
    geom_line(data=data, aes(x=X, y=value, color=variable), size=0.5) +
    geom_abline(slope = 1, intercept = 0, linetype='dashed') +
    scale_color_discrete(name=expression(alpha), labels=paste(x@alpha)) +
    xlab('x') + 
    ylab('Y') 

  p
  return(p)
}

#forec: vector of forecasts, obs: vector of observations
#N.grid.x: length of x-grid to extrapolate over
#gmm: logical, If TRUE: Additionally estimate quantile level parametrically with PointFore
#extrapolate: logical if TRUE: take a grid of x-values between min(forec) and  max(forec)
              #Otherwise take forecast values
estimate_quantile_level <- function(forec, obs, N.grid.x=100, gmm=TRUE, extrapolate=TRUE) {
  #Define grids for alpha and x
  alpha.grid <- seq(0.01, 0.99, 0.01)  ##grid of quantile levels
  
  if (extrapolate==TRUE) {
    x.grid <- seq(min(forec),max(forec),length.out = N.grid.x) ##grid of x values
  }
  
  else{
    x.grid <- forec
  }
  
  #Fit IDR model
  fit <- idr(
    y = obs,
    X = data.frame(X=forec)
  )
  
  #Prediction given forecasts x in x.grid
  prediction <- predict(fit, data = data.frame(X=x.grid))
  
  #Quantiles of Y|x
  #quantiles: matrix of dimension Tx99 
  #rows: x_j, cols: quantiles 0.01,0.02,...,0.99
  quantiles <- qpred(prediction, quantiles = alpha.grid)
  
  alpha_x <- numeric(N.grid.x)
  for (i in 1:length(x.grid)) {
    alpha_x[i] <- alpha.grid[max(which(quantiles[i,]<=x.grid[i]))]
  }
  
  if (gmm==FALSE) {
    res <- new('estimated_levels', x.grid = x.grid, alpha_x=alpha_x, 
               alpha_gmm=NaN, gmm=gmm, res_parametric=NaN, forec=forec)
    return(res)
  }
  
  if (gmm==TRUE){
    res_parametric <- estimate.functional(iden.fct = PointFore::quantiles,
                                          X = forec,
                                          Y = obs,
                                          stateVariable = forec,
                                          instruments = c("lag(Y,2)", "X"),
                                          theta0 = c(0,0),
                                          model = probit_linear)
    
    theta1 <- summary(res_parametric)$coefficients[1,1]
    theta2 <- summary(res_parametric)$coefficients[2,1]
    
    
    alpha_gmm <- probit_linear(x.grid, theta = c(theta1, theta2))
    
    res <- new('estimated_levels', x.grid = x.grid, alpha_x=alpha_x, 
               alpha_gmm=alpha_gmm, gmm=gmm, res_parametric=res_parametric, forec=forec)
    
    return(res)
  }
}


#x object of class estimated_levels
#pdf: logical, if TRUE marginal pdf of x is plotted
#conf_level: vector specifiying confident intervall
#hline: logical, if TRUE hotizontal line at 0.5 is plotted
plot.quantile_level <- function(x, pdf=TRUE, conf.levels=c(0.6,0.9), hline=TRUE) {
  interval_state <- seq(quantile(x@res_parametric$stateVariable, probs = 0.01),quantile(x@res_parametric$stateVariable, probs = 0.99), length.out=100)
  limits <- interval_state[c(1,length(interval_state))]
  
  
  data <- data.frame('x'=x@x.grid, 'alpha_x'=x@alpha_x, 'alpha_gmm'=x@alpha_gmm)
  
  data <- melt(data, id.vars = c('x'))

  p <- ggplot() + 
    geom_line(data = data, aes(x=x, y=value, col=variable), size=0.5) +
    ylim(0,1) +
    xlim(limits) +
    ylab(expression(alpha(x)))

  
  if (x@gmm==FALSE) {
    p <- p + theme(legend.position = "none")
  }
  
  if(pdf==TRUE) {
    p <- p + geom_density(data = data.frame(x@forec),
                          aes(x=x@forec, y=..scaled..),
                          fill = 'green',
                          alpha=.2,
                          show.legend = FALSE) +
      stat_density(aes(x@forec, y=..scaled.., col='pdf estimate of x'),
                       geom='line', position='identity')
  }
  
  if(hline==TRUE) {
    p <- p + geom_hline(yintercept = 0.5, linetype=2)
  }
  
  #conf intervalls
  theta <- c(coef(x@res_parametric$gmm))
  var_theta <- x@res_parametric$gmm$vcov
  
  theta_random <- MASS::mvrnorm(1000,theta,var_theta)
  
  alpha_int <- numeric(length(interval_state))
  alpha_low <- numeric(length(interval_state))
  alpha_high <- numeric(length(interval_state))
  alpha_low2 <- numeric(length(interval_state))
  alpha_high2 <- numeric(length(interval_state))
  
  for ( i in 1:length(interval_state))
  {
    emp.distr <- apply(theta_random, 1,function(theta) x@res_parametric$model(interval_state[i],theta))
    
    alpha_int[i] <- mean(emp.distr)
    alpha_low[i] <- quantile(emp.distr,probs = (1-conf.levels[1])/2)
    alpha_high[i] <- quantile(emp.distr,probs = 1-(1-conf.levels[1])/2)
    alpha_low2[i] <- quantile(emp.distr,probs = (1-conf.levels[2])/2)
    alpha_high2[i] <- quantile(emp.distr,probs = 1-(1-conf.levels[2])/2)
  }
  
  plot_data <- data.frame(cbind(interval_state,alpha_int, alpha_low,alpha_high, alpha_low2,alpha_high2))
  
  p <- p +
    #geom_line(data=plot_data, aes(x=interval_state, y=alpha_int), size=1.2)+
    geom_ribbon(data=plot_data, aes(x=interval_state,ymin=alpha_low,ymax=alpha_high), alpha=0.4)+
    geom_ribbon(data=plot_data, aes(x=interval_state,ymin=alpha_low2,ymax=alpha_high2), alpha=0.2)
  

  
  #add legend
  p <-  p + scale_color_manual(name='', 
                                values=c('alpha_x'='red', 
                                         'alpha_gmm'='blue', 
                                         'pdf estimate of x'='darkgreen'),
                               labels=c('alpha_gmm' = 'estimated quantile level GMM \nwith 60 and 90 percent \nconfidence intervalls',
                                        'alpha_x' = 'estimated quantile level IDR', 
                                        'pdf estimate of x' = 'pdf estimate of x')) 
  
  return(p)
}

###precipitation
#res1 <- estimate_quantiles(precipitation$X, precipitation$Y)
#plot.quantiles(res1)
#res2 <- estimate_quantile_level(precipitation$X, precipitation$Y)
#plot.quantile_level(res2)


###GDP
#res1 <- estimate_quantiles(GDP$forecast, GDP$observation)
#plot.quantiles(res1)
#res2 <- estimate_quantile_level(GDP$forecast, GDP$observation)
#plot.quantile_level(res2)

###ERA5 - Brasilien
#load('Data/precip_era5_025_025_24_2007-2018_-53.5_-6.Rdata')
#res1 <- estimate_quantiles(precip$HRES, precip$OBS*1000)
#plot.quantiles(res1)
#res2 <- estimate_quantile_level(precip$HRES, precip$OBS*1000)
#plot.quantile_level(res2)

###TRMM - Brasilien
load('Data/preD_TRMM_025_025_ECMWF_24_1998-2017_-53.375_-5.875.Rdata')
T <- length(precipData$HRES)
res1 <- estimate_quantiles(precipData$HRES[2:T], precipData$OBS[2:T])
plot.quantiles(res1)
res2 <- estimate_quantile_level(precipData$HRES[2:T], precipData$OBS[2:T])
plot.quantile_level(res2)

