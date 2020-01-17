rm(list=ls())

library(isodistrreg)
library(PointFore)
library(isotone)
library(tidyverse)
library(reshape2)
library(gtable)
library(gridExtra)
library(grid)

y <- precipitation$Y[1:100]
x <- precipitation$X[1:100]

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

p1 <- ggplot() + 
  geom_line(data=plot.data, aes(x=X, y=value, color=variable), 
            size=0.5) +
  geom_abline(slope = 1, intercept = 0, linetype='dashed') +
  scale_color_discrete(name=expression(alpha), 
                       labels=c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9')) +
  guides(color = guide_legend(reverse = TRUE)) +
  ylab('y') + 
  ggtitle('Estimated conditional expectiles - Precipitation London')

##Estimate alpha(x)
alpha.grid <- seq(0.01,0.99,0.01)
N.grid.x <- 100
x.grid <- seq(min(x), max(x), length.out=N.grid.x)
alpha_x <- numeric(N.grid.x)
expectiles <- matrix(ncol=99, nrow=T)

j=1
for(i in alpha.grid) {
  print(i)
  fit.expectiles <- activeSet(total_order, "asyLS", y = y, weights = w, 
                              aw = i, bw = 1-i)
  expectiles[,j] <- fit.expectiles$x
  j=j+1
}

for (i in 1:N.grid.x) {
  alpha_x[i] <- alpha.grid[max(which(expectiles[i,]<=x.grid[i]))]
}

#Estimate expectile level with GMM parametrically
res_parametric <- estimate.functional(iden.fct = PointFore::expectiles,
                                      X = x,
                                      Y = y,
                                      stateVariable = x,
                                      instruments = c("lag(Y,2)", "X"),
                                      theta0 = c(0,0),
                                      model = probit_linear)

theta1 <- summary(res_parametric)$coefficients[1,1]
theta2 <- summary(res_parametric)$coefficients[2,1]
level <- probit_linear(x.grid, theta = c(theta1, theta2))

res <- data.frame('x' = x.grid, 'alpha_x'=alpha_x, 'level'=level)

res <- melt(res, id.vars = c('x'))

p2 <- ggplot() + 
  geom_line(data = res, aes(x=x, y=value, col=variable), size=0.5) +
  #xlim(0,30) +
  ylim(0,1) +
  ylab(expression(alpha(x))) +
  scale_color_discrete(name='Estimated expectile level', labels=c('IDR', 'GMM')) 

#Arrange Plots
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g <- rbind(g1, g2, size="first")
g$widths <- unit.pmax(g1$widths, g2$widths)
grid.newpage()
grid.draw(g)


