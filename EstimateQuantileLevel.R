rm(list=ls())

library(PointFore)
#devtools::install_github("AlexanderHenzi/isodistrreg", build_vignettes = TRUE)
library(isodistrreg)
library(tidyverse)
library(reshape2)

#Formate data sets
dates <-seq(as.Date("2011-07-01"), as.Date("2017/06/30"), "days")
precip <- data.frame('dates'=dates, 'X'=precipitation$X, 'Y'= precipitation$Y)
gdp <- data.frame('dates'=GDP$date, 'X'=GDP$forecast, 'Y'=GDP$observation)

#Choose data set from PointFore 'precip' or 'gdp' (formatted above) and title
#data <- precip
#ti <- 'Estimated quantile level - Precipitation London'
data <- gdp
ti <- 'Estimated quantile level - GDP'

#Initialize variables
T <- length(data$X)  ##data set length

N.grid.x=100
alpha.grid <- seq(0.01, 0.99, 0.01)  ##grid of quantile levels
x.grid <- seq(min(data$X),max(data$X),length.out = N.grid.x) ##grid of forecast levels
x <- data$X ##forecast
alpha_x <- numeric(N.grid.x)  ##vector to save alpha(x)


#Fit IDR model
fit <- idr(
  y = data$Y,
  X = data['X']
)

#Prediction given forecasts x
prediction <- predict(fit, data = data.frame(X=x.grid))


#Quantiles of Y|x
#quantiles: matrix of dimension Tx99 
#rows: x_j, cols: quantiles 0.01,0.02,...,0.99
quantiles <- qpred(prediction, quantiles = alpha.grid)

#for each x in grid find max index where quantile is <= x_i and take the respective alpha
for (i in 1:length(x.grid)) {
  alpha_x[i] <- alpha.grid[max(which(quantiles[i,]<=x.grid[i]))]
}



#Estimate quantile level with GMM parametrically
res_parametric <- estimate.functional(iden.fct = PointFore::quantiles,
                                      X = data$X,
                                      Y = data$Y,
                                      stateVariable = data$X,
                                      instruments = c("lag(Y,2)", "X"),
                                      theta0 = c(0,0),
                                      model = probit_linear)

theta1 <- summary(res_parametric)$coefficients[1,1]
theta2 <- summary(res_parametric)$coefficients[2,1]


level <- probit_linear(x.grid, theta = c(theta1, theta2))

#save result in data frame and reshape for ggplot
res <- data.frame('x' = x.grid, 'alpha_x'=alpha_x, 'level'=level)

res <- melt(res, id.vars = c('x'))

#Plot
ggplot() + 
  geom_line(data = res, aes(x=x, y=value, col=variable), size=0.5) +
  #xlim(0,30) +
  ylim(0,1) +
  ylab(expression(alpha(x))) +
  scale_color_discrete(name='Estimated quantile level', labels=c('IDR', 'GMM')) +
  ggtitle(ti)
