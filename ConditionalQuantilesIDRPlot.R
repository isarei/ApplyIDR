rm(list=ls())

library(PointFore)
library(isodistrreg)
library(tidyverse)
library(reshape2)

#Formate data sets
dates <-seq(as.Date("2011-07-01"), as.Date("2017/06/30"), "days")
precip <- data.frame('dates'=dates, 'X'=precipitation$X, 'Y'= precipitation$Y)
gdp <- data.frame('dates'=GDP$date, 'X'=GDP$forecast, 'Y'=GDP$observation)

#Choose data set from PointFore 'precip' or 'gdp' (formatted above)
#Choose corresponding plot title
data = precip
ti <- "Estimated quantiles of Y|x - Precipitation over London"
#data = gdp
#ti <- "Estimated quantiles of Y|x - GDP data"


#Fit IDR model
fit <- idr(
  y = data$Y,
  X = data['X']
)

#Prediction given forecasts x
prediction <- predict(fit, data = data['X'])

#Quantiles of Y|x
quantiles <- qpred(prediction, quantiles = seq(0.1,0.9,0.1))

#Bring data together and reshape for ggplot2
data <- cbind(data, quantiles)
data <- melt(data, id.vars = c('dates', 'X', 'Y'))

#Plot
p <- ggplot() + 
  geom_line(data=data, aes(x=X, y=value, color=variable), size=0.5) +
  geom_abline(slope = 1, intercept = 0, linetype='dashed') +
  scale_color_discrete(name=expression(alpha), labels=c('0.1', '0.2', '0.3', '0.4', '0.5',
                                                        '0.6', '0.7', '0.8', '0.9')) +
  xlab('x') + 
  ylab('Y') +
  ggtitle(ti)

p