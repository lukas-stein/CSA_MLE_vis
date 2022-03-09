library(miscTools)
library(maxLik)
library(scales)
library(gridExtra)
library(grid)
library(knitr)
library(latex2exp)
source("csa_functions.R")
loc = 0
gam_1 = 1
i=0


while (i <= 3){
  s = round(rnorm(1,0,90000))
  set.seed(s)
  
  csample = rcauchy(10,0,1)
  loc_grid = seq(-750,750,length.out=501)
  ll_value = unlist(lapply(loc_grid, loglike_cauy, gam_val=1, x_vec=csample))
  
  loc_length = 501
  loc_range <- seq(loc-750,loc+750,length.out=loc_length)
  
  csample_ll <- unlist(lapply(loc_range, loglike_cauy, gam_val=gam_1
                              , x_vec=csample))
  minima = c()
  for (i in 2:(length(csample_ll)-1)){
    if (csample_ll[i] <= csample_ll[i-1] && csample_ll[i] <= csample_ll[i+1]){
      minima = c(minima, loc_range[i])
    }
  }
  i = length(minima)
  
  plot(loc_grid, ll_value, type="l", xlab=expression(theta), 
       ylab="-loglike", cex=5, col="green", 
       main=paste("likelihood value (n=7)", sep=""))
}

print(s)

