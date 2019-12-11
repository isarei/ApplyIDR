rm(list=ls())

library(tidyverse)
library(PointFore)
library(reshape2)
library(erboost)

#Data
x <- precipitation$X
y <- precipitation$Y

#Initialize variables
N.grid.x=100
tau.grid <- seq(0.01, 0.99, 0.01)  ##grid of quantile levels
x.grid <- seq(min(x),max(x),length.out = N.grid.x) 
tau_x <- numeric(N.grid.x)

T=length(x)
expectiles <- data.frame(nrow=T)

#Estimate expectiles
for (tau in tau.grid) {
  obj <- erboost(formula = y~x,
                 distribution = list(name="expectile",alpha=tau))
  expectiles <- cbind(expectiles, obj$fit)
} 

expectiles <- expectiles[2:length(tau.grid)]

#solve x=e_{tau(x)}
for (i in 1:length(x.grid)) {
  tau_x[i] <- tau.grid[max(which(expectiles[i,]<=x.grid[i]))]
}

res <- data.frame('x' = x.grid, 'tau_x'=tau_x)
res <- melt(res, id.vars = c('x'))

p <- ggplot() + 
  geom_line(data = res, aes(x=x, y=value, col=variable), size=0.5) +
  #xlim(0,30) +
  ylim(0,1) +
  ylab(expression(tau(x))) +
  scale_color_discrete(name='Estimated expectile level', labels=c('erboost')) 
  #ggtitle(ti)

p
