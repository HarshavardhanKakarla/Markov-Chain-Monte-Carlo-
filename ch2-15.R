set.seed(19)
library(MASS)
library(mvtnorm)

# Metropolis Hastings Algorithm
mh_sampler_mv<-function(start,n_samples=1000,h,dim=100)
{
  samples<-matrix(0,nrow=n_samples,ncol = dim)
  samples[1, ]<-start
  accept<-0
  
  for(i in 2:n_samples)
  {
    proposal <- mvrnorm(1,mu=samples[i-1, ], Sigma=diag(rep(h,dim)))
    
    target_current<-dmvnorm(samples[i-1, ], mean=rep(0,dim), sigma=diag(rep(1,dim)), log = TRUE)
    target_proposal<-dmvnorm(proposal, mean=rep(0,dim), sigma=diag(rep(1,dim)), log=TRUE)
    
    # for numeric stability
    log_alpha <- target_proposal - target_current
    alpha<-min(1, exp(log_alpha))
    
    # Accept-Reject Step
    if(runif(1)<alpha)
    {
      samples[i,]<-proposal
      accept<-accept + 1
    }else
    {
      samples[i,]<-samples[i-1,]
    }
  }
  print(paste("Acceptance rate for h=", h, "is", accept/n_samples))
  return(samples)
}

# For smaller dimension the acceptance rate is non-zero
# For larger dimension the acceptance rate becomes zero, if variance is small we can get non-zero acceptance rate
mh_mv_1<-mh_sampler_mv(start=rep(0,100), n_samples=1000, h=1)
mh_mv_2<-mh_sampler_mv(start=rep(0,100), n_samples=1000, h=5)
mh_mv_3<-mh_sampler_mv(start=rep(0,100), n_samples=1000, h=10)
mh_mv_4<-mh_sampler_mv(start=rep(0,100), n_samples=1000, h=0.01)
mh_mv_5<-mh_sampler_mv(start=rep(0,100), n_samples=1000, h=0.05)

par(mfrow = c(5,1))
# Trace Plots
plot.ts(mh_mv_1[ ,1], main = "Trace Plot for h=1",col = "red", ylab = "Value", xlab = "Iteration")
plot.ts(mh_mv_2[ ,1], main = "Trace Plot for h=5",col = "blue", ylab = "Value", xlab = "Iteration")
plot.ts(mh_mv_3[ ,1], main = "Trace Plot for h=10",col = "orange", ylab = "Value", xlab = "Iteration")
plot.ts(mh_mv_4[ ,1], main = "Trace Plot for h=0.01",col = "red", ylab = "Value", xlab = "Iteration")
plot.ts(mh_mv_5[ ,1], main = "Trace Plot for h=0.05",col = "red", ylab = "Value", xlab = "Iteration")
