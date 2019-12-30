rm(list=ls())

library(PointFore)
library(isodistrreg)
library(tidyverse)
library(reshape2)
library(car)           #for function linearHypothesis (calculate p_wald)
library(grid)          #for adding text to plot with grob
library(cdfquantreg)   #for function scaleTR (transformation of values to unit interval)
library(ggnewscale)    #allows for more than one legend in plot
library(gridExtra)
library(gtable)

setClass("IDR_quantiles", slots=list(X="numeric", Y="numeric", quantiles="matrix", alpha="numeric"))
setOldClass('pointfore')
setClass("estimated_levels", slots=list(x.grid="numeric", alpha_x="numeric", alpha_gmm="numeric", 
                                        gmm="logical", res_parametric='pointfore', forec="numeric",
                                        percentiles="matrix", obs="numeric"))

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
#title: string; plot title
plot.quantiles <- function(x, title) {
  #data
  data <- data.frame('X'=x@X, 'Y'=x@Y)
  data <- cbind(data, x@quantiles)
  data <- melt(data, id.vars = c('X', 'Y'))
  
  #limits
  interval_state <- seq(quantile(x@X, probs = 0.01),
                        quantile(x@X, probs = 0.99), 
                        length.out=100)
  xlimits <- interval_state[c(1,length(interval_state))]

  p <- ggplot() + 
    geom_line(data=data, aes(x=X, y=value, color=variable), 
              size=0.5) +
    geom_abline(slope = 1, intercept = 0, linetype='dashed') +
    scale_color_discrete(name=expression(alpha), 
                         labels=paste(x@alpha)) +
    guides(color = guide_legend(reverse = TRUE)) + 
    theme(legend.justification = "left") + 
    xlim(xlimits) + 
    ylim(min(x@quantiles[,1]),quantile(x@quantiles[,-1],0.99)+5) + 
    xlab('x') + 
    ylab('Y') +
    ggtitle(title)

  p
  return(p)
}

#forec: vector of forecasts, obs: vector of observations
#N.grid.x: length of x-grid to extrapolate over
#gmm: logical, If TRUE: Additionally estimate quantile level parametrically with PointFore
#extrapolate: logical if TRUE: take a grid of x-values between min(forec) and  max(forec)
              #Otherwise take forecast values
#percentiles: liogical if TRUE: Calculate and return conditional quantiles at 0.1,...,0.9
estimate_quantile_level <- function(forec, obs, N.grid.x=100, gmm=TRUE, extrapolate=TRUE, percentiles=TRUE) {
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
    probit0 <- function(stateVariable,theta) probit_linear(stateVariable, theta)*(stateVariable>0)
    res_parametric <- estimate.functional(iden.fct = PointFore::quantiles,
                                          X = forec,
                                          Y = obs,
                                          stateVariable = forec,
                                          instruments = c("lag(Y,2)", "X"),
                                          theta0 = c(0,0),
                                          model = probit0)
                                          #model = probit_linear)
    
    theta1 <- summary(res_parametric)$coefficients[1,1]
    theta2 <- summary(res_parametric)$coefficients[2,1]
    
    
    alpha_gmm <- probit_linear(x.grid, theta = c(theta1, theta2))
    
    if (percentiles==TRUE) {
      prediction2 <- predict(fit, data = data.frame(X=forec))
      perc <- qpred(prediction2, quantiles = seq(0.1,0.9,0.1))
    }
    
    res <- new('estimated_levels', x.grid = x.grid, alpha_x=alpha_x, 
               alpha_gmm=alpha_gmm, gmm=gmm, res_parametric=res_parametric, 
               forec=forec, percentiles=perc, obs=obs)
    
    return(res)
  }
}


#x object of class estimated_levels
#pdf: logical, if TRUE marginal pdf of x is plotted
#conf_level: vector specifiying confident intervall
#hline: logical, if TRUE hotizontal line at 0.5 is plotted
#percentiles: logical, of TRUE add conditional quantile plot
plot.quantile_level <- function(x, pdf=TRUE, conf.levels=c(0.6,0.9), 
                                hline=TRUE, percentiles=TRUE, add.info=TRUE) {
  interval_state <- seq(quantile(x@res_parametric$stateVariable, probs = 0.01),
                        quantile(x@res_parametric$stateVariable, probs = 0.99), 
                        length.out=100) 
  #interval_state <- seq(min(x@forec), max(x@forec), length.out=100)
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
    p <- p + geom_density(data = data.frame(x@forec[x@forec>0]),
                          aes(x=x@forec[x@forec>0], y=..scaled..),
                          fill = 'green',
                          alpha=.2,
                          show.legend = FALSE) +
      stat_density(aes(x@forec[x@forec>0], y=..scaled.., col='pdf estimate of x'),
                       geom='line', position='identity')
  }
  
  if(hline==TRUE) {
    p <- p + geom_hline(yintercept = 0.5, linetype=2)
  }
  
  #Add legend 1
  p <- p + scale_color_manual(name='', 
                       values=c('alpha_x'='red', 
                                'alpha_gmm'='blue', 
                                'pdf estimate of x'='darkgreen'),
                       labels=c('alpha_gmm' = 'estimated quantile level GMM \nwith 60 and 90 percent \nconfidence intervals',
                                'alpha_x' = 'estimated quantile level IDR', 
                                'pdf estimate of x' = 'pdf estimate of x')) +
    theme(legend.justification = "left")
    #theme(panel.background = element_rect(fill = "gray90"),
    #      legend.position=c(0.76,0.3),
    #      legend.position = 'right',
    #      legend.title = element_blank(),
    #      legend.background = element_rect(color = "black", fill = "gray90", 
    #                                       size = 0.5, linetype = "solid"),
    #      legend.key = element_rect(fill = "gray90", color = NaN),
    #      legend.key.width = unit(0.9,"line")) +
    #new_scale('color')
  
  
  if(percentiles==TRUE) {
    data_perc <- data.frame('X'=x@forec, 'Y'=x@obs)
    data_perc <- cbind(data_perc, x@percentiles)
    data_perc <- melt(data_perc, id.vars = c('X', 'Y'))
    p <- p + 
      geom_line(data=data_perc, aes(x=X, y=scaleTR(value), color=variable), 
                size=0.5) +
      #geom_abline(slope = 1, intercept = 0, linetype='dashed') +
      scale_color_discrete(name=expression(alpha), 
                           labels=paste(seq(0.1,0.9,0.1))) +
      guides(color = guide_legend(reverse = TRUE)) + 
      theme(legend.position = 'right') 
      #xlab('x')
      #ylab('Y') 
  }
  
 
  if (add.info==TRUE) {
    #Calculate proportion of zeros in x and y and p-values
    #zeros_X <- round((sum(na.omit(x@forec)==0)/length(x@forec))*100, digits=2)
    zeros_X <- round((sum(na.omit(x@forec)==0)/length(x@forec)), digits=2)
    #zeros_Y <- round((sum(na.omit(x@obs)<10^-5)/length(x@obs))*100, digits=2)
    zeros_Y <- round((sum(na.omit(x@obs)<10^-5)/length(x@obs)), digits=2)
    p_opt <- round(summary(x@res_parametric)$Jtest$test[1,2], digits=2)
    p_wald <- round(unname(attr(linearHypothesis(x@res_parametric$gmm, 'Theta[2]=0'), 'value')[1,1]),
                    digits=2)
    bias <- mean(na.omit(x@forec))-mean(na.omit(x@obs))
    
    
    #Create Text
    grob_p_opt <- grobTree(textGrob(bquote(p[J] * ' = ' * .(p_opt)), x=0.1, y=0.95, hjust=0,
                           gp=gpar(col="black", fontsize=10, fontface="italic")))
    grob_p_wald <- grobTree(textGrob(bquote(p[W] * ' = ' * .(p_wald)), x=0.1, y=0.9, hjust=0,
                                    gp=gpar(col="black", fontsize=10, fontface="italic")))
    #grob_pointmass_x <- grobTree(textGrob(bquote('% 0 in x: ' * .(zeros_X)), x=0.1,  y=0.85, hjust=0,
    #                                      gp=gpar(col="black", fontsize=10, fontface="italic")))
    #grob_pointmass_y <- grobTree(textGrob(bquote('% 0 in y: ' * .(zeros_Y)), x=0.1,  y=0.8, hjust=0,
    #                                      gp=gpar(col="black", fontsize=10, fontface="italic")))
    grob_bias <- grobTree(textGrob(bquote('bias = ' * .(bias)), x=0.1,  y=0.85, hjust=0,
                                   gp=gpar(col="black", fontsize=10, fontface="italic")))
    
    #Add Text to Plot
    p <- p + 
      annotation_custom(grob_p_opt) +
      annotation_custom(grob_p_wald) +
      annotation_custom(grob_bias)
      #annotation_custom(grob_pointmass_x) +
      #annotation_custom(grob_pointmass_y)
    
    #Add Point Mass in 0
    p <- p + 
      geom_point(aes(x=0, y=zeros_X), color='darkgreen', shape=19, size=1) + 
      geom_point(aes(x=0, y=zeros_Y), color='yellow4', shape=19, size=1)
      #geom_segment(aes(x=0, xend=0, y=0, yend=zeros_X), color='darkgreen') +
      #geom_segment(aes(x=0, xend=0, y=0, yend=zeros_Y), color='darkblue')
      
    
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
  #p <-  p + scale_color_manual(name='', 
  #                            values=c('alpha_x'='red', 
  #                                     'alpha_gmm'='blue', 
  #                                     'pdf estimate of x'='darkgreen'),
  #                             labels=c('alpha_gmm' = 'estimated quantile level GMM \nwith 60 and 90 percent \nconfidence intervalls',
  #                                      'alpha_x' = 'estimated quantile level IDR', 
  #                                      'pdf estimate of x' = 'pdf estimate of x')) +
  #  theme(panel.background = element_rect(fill = "gray90"),
  #        legend.position=c(0.76,0.3),
  #        legend.title = element_blank(),
  #        legend.background = element_rect(color = "black", fill = "gray90", 
  #                                         size = 0.5, linetype = "solid"),
  #        legend.key = element_rect(fill = "gray90", color = NaN),
  #        legend.key.width = unit(1.2,"line")) +
  #  new_scale("color") + 
  #  scale_color_discrete(name=expression(alpha), 
  #                       labels=paste(seq(0.1,0.9,0.1))) +
  #  guides(color = guide_legend(reverse = TRUE)) + 
    
  
  return(p)
}




###precipitation
#res1 <- estimate_quantiles(precipitation$X, precipitation$Y)
#p1 <- plot.quantiles(res1,
#                     title='Estimated conditional quantiles and quantile level \nLocation: London')
#res2 <- estimate_quantile_level(precipitation$X, precipitation$Y)
#p2 <- plot.quantile_level(res2, percentiles = FALSE)


###GDP
#res1 <- estimate_quantiles(GDP$forecast, GDP$observation)
#plot.quantiles(res1)
#res2 <- estimate_quantile_level(GDP$forecast, GDP$observation)
#plot.quantile_level(res2)

###ERA5
##Brasilien
#load('Data/precip_era5_025_025_24_2007-2018_-53.5_-6.Rdata')
##Kongo
load('Data/precip_era5_025_025_24_2007-2018_24_-7.5.Rdata')

res1 <- estimate_quantiles(precip$HRES, precip$OBS*1000)
p1 <- plot.quantiles(res1, 
                     title='Estimated conditional quantiles and quantile level \nLocation: Kongo (24,-7.5)')
res2 <- estimate_quantile_level(precip$HRES, precip$OBS*1000)
p2 <- plot.quantile_level(res2, percentiles = FALSE)


###TRMM - Brasilien
#load('Data/preD_TRMM_025_025_ECMWF_24_1998-2017_-53.375_-5.875.Rdata')
#T <- length(precipData$HRES)
#res1 <- estimate_quantiles(precipData$HRES[2:T], precipData$OBS[2:T])
#plot.quantiles(res1)
#res2 <- estimate_quantile_level(precipData$HRES[2:T], precipData$OBS[2:T])
#plot.quantile_level(res2)

#Arrange Plots
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g <- rbind(g1, g2, size="first")
g$widths <- unit.pmax(g1$widths, g2$widths)
grid.newpage()
grid.draw(g)




