" Code: Consider target distribution F = N(0,1) and proposal Q(x,·) = N(x,h) (this
is obviously not a realistic scenario). Run a MH algorithm with starting value X0 = 0
and obtain 1000 samples for each of h = 1,5,10. What are the estimated acceptance
probabilities for each h? Which h seems best from the trace plots? "

set.seed(19)
library(mcmcse)

# for plotting
foo.x = seq(-5, 5, length = 100)
foo.y = dnorm(seq(-5, 5, length = 100))

# MH sampler
mh_sampler<-function(start=0,n,h=1)
{
  x <- numeric(length = n)
  x[1]<-1
  accept<-0
  for(i in 2:n)
  {
    proposal<-rnorm(1, mean = x[i-1], sd=sqrt(h))
    
    target_current <- dnorm(x[i-1], mean = 0, sd = 1)
    target_proposal <- dnorm(proposal, mean = 0, sd = 1)
    # Symmetric Distribution
    alpha <- min(1, target_proposal / target_current)
    
    if(runif(1)<alpha)
    {
      ## Accept step
      x[i]<-proposal
      accept<-accept+1
    }
    else{
      ## Rejection Step
      x[i]<-x[i-1]
    }
  }
    print(paste("Acceptance rate is", accept/n))
    return(x)
}

mh1 <- mh_sampler(start=0, n=1000, h=1)
mh2 <- mh_sampler(start=0, n=1000, h=5)
mh3 <- mh_sampler(start=0, n=1000, h=10)

par(mfrow = c(3,3))

# Trace Plots
plot.ts(mh1, type="l", col="red", main="Trace Plot for h=1",  xlab="Iteration", ylab="X")
plot.ts(mh2, type="l", col="blue", main="Trace Plot for h=5",  xlab="Iteration", ylab="X")
plot.ts(mh3, type="l", col="orange", main="Trace Plot for h=10",  xlab="Iteration", ylab="X")

# ACF Plots
acf(mh1, main="h=1")
acf(mh2,main="h=5")
acf(mh3,main="h=10")

# Density Plots
plot(density(mh1), col = "red", main = "h = 1")
lines(foo.x, y = foo.y, col = "black")
plot(density(mh2), col = "blue", main = "h = 5")
lines(foo.x, y = foo.y, col = "black")
plot(density(mh3), col = "blue", main = "h = 10")
lines(foo.x, y = foo.y, col = "black")