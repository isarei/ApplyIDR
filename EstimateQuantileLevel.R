rm(list=ls())

library(PointFore)
#devtools::install_github("AlexanderHenzi/isodistrreg", build_vignettes = TRUE)
library(isodistrreg)
library(tidyverse)
library(reshape2)

dates <-seq(as.Date("2011-07-01"), as.Date("2017/06/30"), "days")
precip <- data.frame('dates'=dates, 'X'=precipitation$X, 'Y'= precipitation$Y)

#Choose data set from PointFore 'precip' or 'gdp' (formatted above) and title
data <- precip
ti <- 'Estimated quantile level - Precipitation London'
#data <- gdp
#ti <- 'Estimated quantile level - GDP'

#Initialize variables
T <- length(data$X)  ##data set length
alphas <- seq(0.01, 0.99, 0.01)  ##quantile levels
x <- data$X ##forecast
alpha_x <- numeric(T)  ##vector to save alpha(x)

#Fit IDR model
fit <- idr(
  y = data$Y,
  X = data['X']
)

#Prediction given forecasts x
prediction <- predict(fit, data = data['X'])

#Quantiles of Y|x
#quantiles: matrix of dimension Tx99 
#rows: x_j, cols: quantiles 0.01,0.02,...,0.99
quantiles <- qpred(prediction, quantiles = alphas)

#for each x_i (row) find max index where quantile is <= x_i and take the respective alpha
for (i in seq(1, T, 1)) {
  alpha_x[i] <- alphas[max(which(quantiles[i,]<=x[i]))]
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

level <- probit_linear(res_parametric$stateVariable, theta = c(theta1, theta2))

#save result in data frame and reshape for ggplot
res <- data.frame('x' = x[3:T], 'alpha_x'=alpha_x[3:T], 'level'=level)
res <- melt(res, id.vars = c('x'))

#Plot
ggplot() + 
  geom_line(data = res, aes(x=x, y=value, col=variable), size=0.5) +
  #xlim(0,30) +
  ylim(0,1) +
  ylab(expression(alpha(x))) +
  scale_color_discrete(name='Estimated quantile level', labels=c('IDR', 'GMM')) +
  ggtitle(ti)

