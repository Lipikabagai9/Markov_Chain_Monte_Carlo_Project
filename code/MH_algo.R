# Metropolis-Hastings for Mixture of Gaussian Target Distribution 
target_density = function(x)
{
  mu1 = -2; mu2 = -1; mu3 = 1; mu4 = 2
  sigma1 = 0.05; sigma2 = 0.25; sigma3 = 0.25; sigma4 = 0.1
  density_val = 0.2 * dnorm(x, mean=mu1, sd=sigma1) + 
    0.2 * dnorm(x, mean=mu2, sd=sigma2) + 
    0.3 * dnorm(x, mean=mu3, sd=sigma3) + 
    0.3 * dnorm(x, mean=mu4, sd=sigma4) 
  return(density_val)
}

# 1. Uniform density proposal
metropolis_hastings = function(n, initial_state, step_size){
  samples = numeric(n)
  samples[1] = initial_state
  
  current_x = initial_state
  accepted_count = 0
  
  for (i in 2:n){
    # Using uniform proposal
    noise = runif(1, min = -step_size, max = step_size)
    proposed_x = current_x + noise
    
    # Calculate Acceptance Probability (alpha)
    # Proposal is symmetric, so q(x|y) cancels out
    ratio = target_density(proposed_x) / target_density(current_x)
    alpha = min(1, ratio)
    
    if(runif(1) < alpha){
      current_x = proposed_x
      accepted_count = accepted_count + 1
    } 
    samples[i] = current_x
  }
  
  # Acceptance rate
  acc_rate = accepted_count / (n - 1)

  return(list(chain = samples, rate = acc_rate))
}


##########################################
# 2. Using Normal Distribution proposal centered at current_x
metropolis_hastings_gaussian = function(n, initial_state, step_size){
  samples = numeric(n)
  samples[1] = initial_state
  
  current_x =initial_state
  accepted_count = 0
  
  for (i in 2:n) {
    proposed_x = rnorm(1, mean = current_x, sd = step_size)
    
    ratio = target_density(proposed_x) / target_density(current_x)
    alpha = min(1, ratio)
    
    if (runif(1) < alpha) {
      current_x = proposed_x     
      accepted_count = accepted_count + 1
    }
    samples[i] <- current_x
  }
  
  acc_rate <- accepted_count / (n- 1)
  return(list(chain = samples, rate = acc_rate))
}

##########################################
# Visualization


n_iter = 10000
par(mfrow = c(1, 2))
# Analysis for different step sizes
step_sizes = c(0.1, 2.0, 50.0) 
titles = c("Small Step (Slow)", "Optimal Step (Good)", "Large Step (Stuck)")

# Uniform proposal
for (i in 1:3) {
  res = metropolis_hastings(n = n_iter, initial_state = 0, step_size = step_sizes[i])
  plot(res$chain, type = "l", col = "darkgreen", 
       main = paste(titles[i], "\nRate:", round(res$rate, 2)),
       ylab = "x", xlab = "Iteration", ylim = c(-3, 3))
  hist(res$chain, probability = TRUE, col = "lightblue", 
       breaks = 40, main = "Estimated vs True Density", xlab = "x")
  curve(target_density(x), col = "red", lwd = 2, add = TRUE)
}

# Gaussian Proposal
for (i in 1:3) {
  res = metropolis_hastings(n = n_iter, initial_state = 0, step_size = step_sizes[i])
  plot(res$chain, type = "l", col = "darkgreen", 
       main = paste(titles[i], "\nRate:", round(res$rate, 2)),
       ylab = "x", xlab = "Iteration", ylim = c(-3, 3))
  hist(res$chain, probability = TRUE, col = "lightblue", 
       breaks = 40, main = "Estimated vs True Density", xlab = "x")
  curve(target_density(x), col = "red", lwd = 2, add = TRUE)
}

############################################

