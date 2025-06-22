set.seed(19)
library(mcmcse)


log_f <- function(x) {
  return(-x^4)
}
grad_log_f <- function(x){
  return(-4*x^3)
}

log.mala.den <- function(curr, goingto, h)
{
  dnorm(goingto, mean = curr + h/2*grad_log_f(curr), sd = sqrt(h), log = TRUE)
}

mala_1d <- function(start, n_iter, step_size) {
  x <- numeric(n_iter)
  x[1] <- start
  accept <- 0
  
  for (i in 2:n_iter) {
    grad <- grad_log_f(x[i-1])
    proposal <- x[i-1] + 0.5*step_size*grad + rnorm(1, 0, sqrt(step_size))
    
    grad_prop <- grad_log_f(proposal)
    log_alpha <- log_f(proposal) - log_f(x[i-1]) +
      log.mala.den(proposal, x[i-1], step_size) -
      log.mala.den(x[i-1], proposal, step_size)
    
    if (log(runif(1)) < log_alpha) {
      x[i] <- proposal
      accept <- accept + 1
    } else {
      x[i] <- x[i-1]
    }
  }
  cat("Start:",start, "Step Size:",step_size, "Acceptance rate:", accept/n_iter,"\n")
  return(list(samples=x, acceptance=accept/n_iter))
}

start_values <- c(0, 10, 100)  #  Starting points
step_sizes <- c(0.1, 0.5, 1.0)  # Different step sizes

results <- list()
# Test with different starting values and step sizes
for (start in start_values) {
  for (step in step_sizes) {
    result <- mala_1d(start, n_iter = 10000, step_size=step)
    results[[paste("Start:", start, "Step:", step)]] <- result
  }
}

par(mfrow = c(3, 3))
for (name in names(results)) {
  plot.ts(results[[name]]$samples, main = name, ylab = "Value")
}

# Observations: 
# - Chain mixes well from reasonable starts
# - From extreme starts (e.g., 10), initial transient period longer
# - Requires careful step_size tuning for good acceptance

#par(mfrow = c(3,1), mar = c(2, 4, 2, 1))
#plot.ts(out1, main="MALA for strating with 0", col="red")
#plot.ts(out2, main="MALA for strating with 0", col="blue")
#plot.ts(out3, main="MALA for strating with 0", col="orange")
