set.seed(19)
library(mcmcse)

## Unadjusted Langevin Algorithm:
# F = N(5, 20) true distribution
# Discretized Langevin dynamics: x_{k+1} = x_k + h*(5 - x_k)/20 + sqrt(2h)*Z
h_values = c(.001,.01,.1,1,10)

ula_algo <- function( h=1, n_iter=10000)
{
  x<-numeric(n_iter)
  x[1]<-0
  
  for(i in 2:n_iter)
  {
    drift <- h*(5-x[i-1])/20
    noise <- sqrt(2*h)*rnorm(1)
    x[i]<-x[i-1] + drift + noise
  }
  x
}

ula_results <- lapply(h_values, ula_algo)

par(mfrow=c(5,1))
for(i in 1:5)
{
  plot(density(ula_results[[i]]), main=paste("ULA h=", h_values[i]))
  curve(dnorm(x,5,sqrt(20)), add=TRUE, col="red")
}

