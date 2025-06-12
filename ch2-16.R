set.seed(19)
library(MASS)
library(mvtnorm)

# MH Algorithm
mh_algorithm<- function(h, n_samples = 1000, dim=5) 
{
  x <- matrix(0, nrow = n_samples, ncol = 5)  
  accept <- 0
  
  for (i in 2:n_samples) 
  {
    proposal <- mvrnorm(1, mu = x[i - 1, ], Sigma = diag(rep(h, dim)))
    
    target_current<-dmvnorm(x[i-1, ], mean=rep(0,dim), sigma = diag(c(1, 2, 5, 10, 100)), log = TRUE)
    target_proposal<-dmvnorm(proposal, mean=rep(0,dim), sigma = diag(c(1, 2, 5, 10, 100)), log = TRUE)
    
    
    # for numeric stability
    log_alpha <- target_proposal - target_current
    acceptance_ratio<-min(1, exp(log_alpha))
    
    if (runif(1) < acceptance_ratio) {
      x[i, ] <- proposal
      accept <- accept + 1
    } else {
      x[i, ] <- x[i - 1, ]
    }
  }
  
  list(samples = x, acceptance_rate = accept / n_samples)
}

# Tuning h for having ~30% acceptance rate
h_values <- seq(0.1, 5, by = 0.1)
optimal_h_value <- 0
best_acceptance <- 0

# Iterating for optimal h value in seq(0.1,5) with step size 0.1
for (h in h_values) 
{
  res <- mh_algorithm(h)
  acc_rate <- res$acceptance_rate
  
  # Printing h acceptance value for each h value
  cat(sprintf("For h=%.3f -> Acceptance rate=%.3f\n", h, acc_rate))
  
  # Condition for checking optimal h value and acceptance value
  if (abs(acc_rate - 0.3) < abs(best_acceptance - 0.3)) {
    optimal_h_value <- h
    best_acceptance <- acc_rate
  }
}

cat("Optimal value of h is :", optimal_h_value, "\n Estimated acceptance rate is:", best_acceptance, "\n")
