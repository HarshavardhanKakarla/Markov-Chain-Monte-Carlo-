set.seed(19)
library(mcmcse)
data(Pima.tr)

X <- model.matrix(~ npreg + glu + bp + skin + bmi + ped + age, data=Pima.tr) 
y <- ifelse(Pima.tr$type == "Yes", 1, 0)

log_post <- function(beta) 
{
  linear_pred <- X %*% beta
  sum(y * linear_pred - log(1 + exp(linear_pred)))
}

grad_log_post <- function(beta)
{
  p <- 1 / (1+exp(-X %*% beta))
  t(X) %*% (y-p)
}

log_mala_den <- function(from, to, step_size, grad_log_fun) {
  # grad_log_fun should be a function that returns gradient at a given point
  mean_cal <- from + 0.5 * step_size * grad_log_fun(from)
  dnorm(to, mean = mean_cal, sd = sqrt(step_size), log = TRUE)
}

mala_bay_log<-function(iter=10000, step_size=0.1)
{
  p <- ncol(X)
  beta <- rep(0,p)
  chain <- matrix(NA, iter, p)
  accept<-0
  
  for(i in 1:iter)
  {
    grad <- grad_log_post(beta)
    proposal <- beta + 0.5*step_size*grad + rnorm(p, 0, sqrt(step_size))
    
    grad_prop<- grad_log_post(proposal)
    
    log_alpha <- log_post(proposal) - log_post(beta) + 
      # used sum function as dnorm gives vectors of log densities to get scalars
       sum(log_mala_den(beta,proposal, step_size,grad_log_post))-
       sum(log_mala_den(proposal, beta, step_size, grad_log_post))
    
    if(runif(1) < log_alpha)
    {
      beta <- proposal
      accept <- accept+1
    }else{
      chain[i,] <- beta
    }
  }
  list(samples=chain, acceptance=accept/iter)
}

step_size <- 0.01
for (i in 1:20) {
  out <- mala_bay_log(2000, step_size)
  cat("Step:", step_size, "Acceptance rate:", out$acceptance, "\n")
  if (out$acceptance < 0.55) step_size <- step_size * 0.8
  else if (out$acceptance > 0.65) step_size <- step_size * 1.2
  else break
}

result_mala <- mala_bay_log(10000, step_size)
