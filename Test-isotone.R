rm(list=ls())

library(PointFore)
library(isotone)
library(tidyverse)
library(reshape2)

x <- precipitation$X[1:100]
y <- precipitation$Y[1:100]

T <- length(x)

total_order <- cbind(1:(T-1), 2:T)

w <- rep(1, T)

##ActiveSet
fit.qua_active <- activeSet(total_order, "quantile", y = y, weights = w, 
                     aw = 0.5, bw = 0.5)
##PAVA
fit.qua_pava <- gpava(z=x, y=y, solver=weighted.fractile, p=0.5)

##IDR
fit <- idr(
  y = y,
  X = data.frame(X=x)
)
prediction <- predict(fit, data = data.frame(X=x))
quantiles <- qpred(prediction, quantiles = 0.5)


res <- data.frame('x'=x, 'activeSet'=fit.qua_active$x, 'pava'=fit.qua_pava$x, 'IDR'=quantiles)
res <- melt(res, id.vars = c('x'))
ggplot() + 
  geom_line(data = res, aes(x=x, y=value, col=variable), size=0.5) +
  scale_color_discrete(name='') + 
  ylab('Y')


