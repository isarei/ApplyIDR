rm(list=ls())

library(PointFore)
#devtools::install_github("AlexanderHenzi/isodistrreg", build_vignettes = TRUE)
library(isodistrreg)
library(tidyverse)

dates <-seq(as.Date("2011-07-01"), as.Date("2017/06/30"), "days")
precip <- data.frame('dates'=dates, 'X'=precipitation$X, 'Y'= precipitation$Y)

#all
T = length(precip$dates)
#for test
#T = 30


fit <- idr(y=precip$Y, X=precip['X'])
alpha_t <- numeric(T)
diff_all <- numeric(T)
for (t in seq(1,T,1)) {
  print(t)
  diff <- 10^6
  for (alpha in seq(0.01,0.99,0.01)) {
    print(alpha)
    
    x <- (precip$X)[t]
    
    pred <- predict(fit, data=data.frame(X=x))
    q <- qpred(pred, quantiles = alpha)
               
    diff_new <- abs(q-x)
    
    if (diff_new < diff){
      diff <- diff_new
      alpha_t[t] <- alpha
    }
  }
  diff_all[t] = diff
}
