set.seed(19)
library(MASS)
library(mvtnorm)

# Target covariance matrix (diagonal: 1, 2, 5, 10, 100)
Sigma_diag <- c(1, 2, 5, 10, 100)
dim <- length(Sigma_diag)

# adaptive MH sampler
mh_adaptive <- function(n_samples = 10000, target_accept = 0.3) 
{
  # Initialize
  current <- rep(0, dim)
  samples <- matrix(0, nrow = n_samples, ncol = dim)
  samples[1, ] <- current
  accept <- 0
  
  # Initialize proposal parameters
  h <- Sigma_diag  
  s <- 2.4^2 / dim  # Optimal scaling factor for multivariate normal
  
  # Main sampling loop
  for (i in 2:n_samples) {
    # Propose new point
    proposal <- current + rnorm(dim, mean = 0, sd = sqrt(s * h))
    
    # Calculate log-target densities
    log_target_current <- dmvnorm(current, mean=rep(0,dim), sigma = diag(c(1, 2, 5, 10, 100)), log = TRUE)
    log_target_proposal <- dmvnorm(proposal, mean=rep(0,dim), sigma = diag(c(1, 2, 5, 10, 100)), log = TRUE)
    
    # Acceptance probability
    log_alpha <- log_target_proposal - log_target_current
    if (log(runif(1)) < log_alpha) {
      current <- proposal
      accept <- accept + 1
    }
    samples[i, ] <- current
    
   
    if (i > 1000 && i %% 500 == 0) {
      # Update empirical variances
      emp_var <- apply(samples[1:i, ], 2, var)
      
      # Update proposal variances (h) and scaling factor (s)
      h <- emp_var
      current_accept <- accept / i
      if(current_accept < target_accept)
      {
        s <- s*0.95
      }else
      {
        s<-s*1.05
      }
    }
  }
  
  # Discard first 1000 samples
  final_samples <- samples[1001:n_samples, ]
  
  list(samples = final_samples,acceptance_rate = accept / (n_samples - 1),tuned_h = s * h)
}

# Run the sampler with corrected function name
result <- mh_adaptive()

# Results
cat("Acceptance rate:", round(result$acceptance_rate, 3), "\n")
cat("Tuned proposal variances \n h1   h2   h3   h4   h5 \n", round(result$tuned_h, 2), "\n")
