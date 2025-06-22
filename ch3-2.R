set.seed(19)
library(MASS)
library(mcmcse)
data(Pima.tr)

X <- model.matrix(~ npreg + glu + bp + skin + bmi + ped + age, data=Pima.tr) # Add intercept
y <- ifelse(Pima.tr$type == "Yes", 1, 0)
head(X)
head(y)
dim(X)

log_posterior <- function(beta) 
{
  linear_pred <- X %*% beta
  sum(y * linear_pred - log(1 + exp(linear_pred)))
}

rwMh <- function(n_iter = 10000, proposal_sd=0.1)
{
  p<-ncol(X)
  current_beta<-rep(0,p)
  chain <- matrix(NA, n_iter, p)
  accept<-0
  for(i in 1:n_iter)
  {
    proposed_beta <- current_beta + rnorm(p, 0, proposal_sd)
    
    log_posterior_current <- log_posterior(current_beta)
    log_posterior_proposed <- log_posterior(proposed_beta)
    
    acceptance_prob <- exp(log_posterior_proposed - log_posterior_current)
    if(runif(1) < acceptance_prob)
    {
      current_beta<-proposed_beta
      accept <- accept+1
    }else {
      chain[i,]<-current_beta
    }
  }
  list(chain=chain, acceptance=accept/n_iter)
  
}

proposal_sd <- 0.1
for (i in 1:20) {
  out <- rwMh(2000, proposal_sd)
  cat("SD:", proposal_sd, "Acceptance rate:", out$acceptance, "\n")
  if (out$acceptance < 0.25) proposal_sd <- proposal_sd * 0.9
  else if (out$acceptance > 0.35) proposal_sd <- proposal_sd * 1.1
  else break
}


