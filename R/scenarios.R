library(tidyverse)

n=500
# Set up variables, x.d = c(x,c) ----
x <- c(rep(0, n/2), rep(1, n/2))
c1 <- c(rep(0, 0.6*n/2), rep(1, 0.4*n/2))
c2 <- c(rep(0, 0.6*n/2), rep(1, 0.4*n/2))
c <- as.matrix(data.frame( age = c(c1, c2))) 
x.d <- cbind(x, c)

mdist = "ZINB"

b0 = -2.0;  bx = 0.30;  bc = 0.2; bm = 0.02; bi = 0.20; bxi = 0; bxm = 0
g0 = 2.1;   gx = -1.5;  gc = 0.4
a0 = 3.0;   ax = 0.25;  ac = 0.35; theta = 5

# data = source("datasim.R")
