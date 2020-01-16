rm(list=ls())

library(isodistrreg)
library(PointFore)
library(isotone)
library(tidyverse)
library(reshape2)

y <- precipitation$Y[1:200]
x <- precipitation$X[1:200]

#P: Sort data such that appropriate to the defined total_order
y <- y[order(x)]
x <- x[order(x)]

#Take length of forecast vector
T <- length(x)

#Define total order
total_order <- cbind(1:(T-1), 2:T)

#take unit weigths
w <- rep(1, T)

#initialize result matrix
expectiles <- matrix(ncol=9, nrow=T)

for(i in seq(0.1,0.9,0.1)) {
  fit.expectiles <- activeSet(total_order, "asyLS", y = y, weights = w, 
                              aw = i, bw = 1-i)
  expectiles[,i*10] <- fit.expectiles$x
}

plot.data <- data.frame('X'=x, 'Y'=y)
plot.data <- cbind(plot.data, expectiles)
plot.data <- melt(plot.data, id.vars = c('X', 'Y'))

p <- ggplot() + 
  geom_line(data=plot.data, aes(x=X, y=value, color=variable), 
            size=0.5) +
  geom_abline(slope = 1, intercept = 0, linetype='dashed') +
  scale_color_discrete(name=expression(alpha), 
                       labels=c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9')) +
  guides(color = guide_legend(reverse = TRUE)) +
  ylab('y') + 
  ggtitle('Estimated conditional expectiles - Precipitation London')
p




