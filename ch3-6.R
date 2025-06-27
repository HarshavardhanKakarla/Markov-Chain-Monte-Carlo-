set.seed(19)
library(mcmcse)

# True distribution: N(5, 20)
h_values <- c(0.001, 0.01, 0.1, 1, 10)
# Gradient of log-target for standard normal N(5, 20)
grad_log_f <- function(x) {
  return((5 - x) / 20)
}

# Generalized log-density of MALA Gaussian proposal
log_mala_den <- function(from, to, h) {
  mean_shift <- from + 0.5 * h * grad_log_f(from)
  dnorm(to, mean = mean_shift, sd = sqrt(2 * h), log = TRUE)
}
# MALA implementation
mala_normal <- function(h, n_iter=10000) {
  x <- numeric(n_iter)
  x[1] <- 0
  accept <- 0
  
  for (i in 2:n_iter) {
    # Proposal
    drift <- h*(5 - x[i-1])/20
    proposal <- x[i-1] + drift + sqrt(2*h)*rnorm(1)
    
    # Reverse move
    drift_rev <- h*(5 - proposal)/20
    log_alpha <- dnorm(proposal, 5, sqrt(20), log=TRUE) - 
      dnorm(x[i-1], 5, sqrt(20), log=TRUE) +
      log_mala_den(x[i-1], proposal+drift_rev,h)-
      log_mala_den(proposal, x[i-1]+drift,h)
    
    if (log(runif(1)) < log_alpha) {
      x[i] <- proposal
      accept <- accept + 1
    } else {
      x[i] <- x[i-1]
    }
  }
  list(samples=x, acceptance=accept/n_iter)
}
mala_results <- lapply(h_values, function(h) mala_normal(h)$samples)

# Plot density estimates
par(mfrow=c(5,1),mar = c(2, 4, 2, 1))
for (i in 1:5) {
  plot(density(mala_results[[i]]), main=paste("MALA h=", h_values[i]))
  curve(dnorm(x,5,sqrt(20)), add=TRUE, col="red")
}




