rm(list=ls())

library(tidyverse)
library(PointFore)
library(reshape2)
library(erboost)

x <- precipitation$X
y <- precipitation$Y

T=length(x)

expectiles <- data.frame(nrow=T)

for (i in seq(0.1,0.9,0.1)) {
  obj <- erboost(formula = y~x,
                distribution = list(name="expectile",alpha=i))
  expectiles <- cbind(expectiles, obj$fit)
}  

data <- data.frame('X'=x, 'Y'=y)
data <- cbind(data, expectiles[,2:10])
data <- melt(data, id.vars = c('X', 'Y'))
p <- ggplot() + 
  geom_line(data=data, aes(x=X, y=value, color=variable), size=0.5) +
  geom_abline(slope = 1, intercept = 0, linetype='dashed') +
  scale_color_discrete(name=expression(tau), labels=paste(seq(0.1,0.9,0.1))) +
  xlab('x') + 
  ylab('Y') 

p
