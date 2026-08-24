# Define Potential
potential = function(x)
{
  mu1 = -2; mu2 = -1; mu3 = 1; mu4 = 2
  sigma1 = 0.05; sigma2 = 0.25; sigma3 = 0.25; sigma4 = 0.1
  density_val = 0.2 * dnorm(x, mean=mu1, sd=sigma1) + 
    0.2 * dnorm(x, mean=mu2, sd=sigma2) + 
    0.3 * dnorm(x, mean=mu3, sd=sigma3) + 
    0.3 * dnorm(x, mean=mu4, sd=sigma4) 
  return(-log(density_val))
}

# define a grid for calculations
grid = seq(-4, 4, length.out = 2000)

# Pre-calculate the Potential on this grid (black line data)
potential_on_grid = sapply(grid, potential)

# Moreau Envelope and Proximal Operator Calculation
# we solve the optimization problem: min_y { U(y) + (1/2t) * ||x - y||^2 }
moreau_func = function(x, t)
{
  if(t==0){
    return(list(val = potential(x), 0))
  }
  penalty = (1/(2*t))*(x-grid)^2
  objectives = potential_on_grid + penalty
  min_idx = which.min(objectives)
  min_value = min(objectives)                
  return(list(val = min_value, prox = grid[min_idx])) # Return the minimized value and proximal  
}

density_val = function(x, t){
  return(exp(-1*moreau_func(x, t)$val))
}


# Parallel Tempering Algorithm
# 1. Using MH-ALgo based Markov chains 
run_parallel_mh = function(iter, n_chains) 
{
  # T=1 is target. T=High is flat.
  temps = exp(seq(log(1), log(50), length.out = n_chains))
  betas = 1/temps 
  
  chain_states = matrix(0, nrow = iter, ncol = n_chains)
  
  #Initialize
  current_x = rep(0, n_chains)
  
  tot_iter = 0
  swap_success = 0
  
  for(i in 1:iter){
    for(k in 1:n_chains){
      # Propose new spot
      x_old = current_x[k]
      x_new = x_old + rnorm(1, mean=0, sd=0.5) # Random walk
      
      # Calculate Energy Difference
      U_old = potential(x_old)
      U_new = potential(x_new)
      
      # Acceptance Ratio = exp(-beta *(U_new - U_old))
      log_ratio = -betas[k] *(U_new - U_old)
      
      if(log(runif(1)) < log_ratio){
        current_x[k] = x_new
      }
    }
    
    # swapping chain k and k+1
    for (k in 1:(n_chains - 1)) {
      chain_i = k
      chain_j = k + 1
      
      x_i = current_x[chain_i]
      x_j = current_x[chain_j]
      
      # The Swap Acceptance Formula = min(1, exp((Beta_i-Beta_j)*(U(x_i)-U(x_j))))
      d_beta = betas[chain_i]-betas[chain_j]
      d_energy = potential(x_i)-potential(x_j)
      
      log_swap_ratio = d_beta * d_energy
      
      if (log(runif(1)) < log_swap_ratio) {
        # Perform Swap
        temp = current_x[chain_i]
        current_x[chain_i] = current_x[chain_j]
        current_x[chain_j] = temp
        
        swap_success = swap_success + 1
      }
      tot_iter = tot_iter + 1
    }
    # Store states
    chain_states[i, ] = current_x
  }
  
  cat("Swap Acceptance Rate:", round(swap_success/tot_iter, 3), "\n")
  return(list(chains = chain_states, temps = temps))
}

# Run and Visualize
res1 = run_parallel_mh(iter = 5000, n_chains = 5)

# Plotting
par(mfrow=c(1,2))

# Trace of the Target Chain
plot(res1$chains[, 1], type="l", col="blue",
     main="Chain 1 (Target Temp T=1)", 
     ylab="x", xlab="Iteration", ylim =c(-3,3))
abline(h=c(-2, -1, 1, 2), col="red", lty=2)


# Plot the Histogram 
hist(res1$chains[,1], probability = TRUE, breaks = 50, col = "lightblue",
     main = "Final Distribution (Parallel MH)",
     ylim = c(0,1.5), xlab = "Position (x)", ylab = "Density")

# Overlay the True Target Density (Red Line)
curve(exp(-potential(x)), from = -4, to = 4, col = "red", lwd = 2, add = TRUE)

# Add legend
legend("topright", legend=c("PT Samples", "True Target"), 
       col=c("lightblue", "red"), lwd=c(10, 2))




# 2. Using MYPT with RW 
run_PT_RW= function(iter, n_chains) {
  # T=0 is target. T=High is flat.
  t_ladder = c(0, exp(seq(log(0.01), log(5), length.out = n_chains-1)))
  
  # Step size
  tau = 0.05
  
  chain_states = matrix(0, nrow = iter, ncol = n_chains)
  current_x = rep(0, n_chains) # Start all at 0
  
  swap_success = 0
  tot_iter = 0
  
  for(i in 1:iter){
    for(k in 1:n_chains){
      t_val = t_ladder[k]
      sigma = 0.5
      
      # RW proposal (symmetric) 
      x_old = current_x[k]
      x_new = x_old + rnorm(1, mean=0, sd=sigma)
      
      # RW accept-reject step
      # Calculate Energy Difference
      U_old = moreau_func(x_old, t_val)$val
      U_new = moreau_func(x_new, t_val)$val
      
      # Acceptance Ratio
      log_ratio = (U_old - U_new)
      
      if(log(runif(1)) < log_ratio){
        current_x[k] = x_new
      }
    }
    
    # Performing the swaps
    for (k in 1:(n_chains - 1)) {
      idx1 = k
      idx2 = k + 1
      
      x1 = current_x[idx1]
      x2 = current_x[idx2]
      t1 = t_ladder[idx1]
      t2 = t_ladder[idx2]
      
      # Ratio = exp(-[(M^t1(x2) + M^t2(x1)) - (M^t1(x1) + M^t2(x2))])
      
      # Calculate Energy of x2 if it moves to chain 1 (and vice versa)
      E_x2_at_t1 = moreau_func(x2, t1)$val
      E_x1_at_t2 = moreau_func(x1, t2)$val
      
      # Calculate Current Energies
      E_x1_at_t1 = moreau_func(x1, t1)$val
      E_x2_at_t2 = moreau_func(x2, t2)$val
      
      # Proposed Energy - Current Energy
      delta_energy = (E_x2_at_t1 + E_x1_at_t2) - (E_x1_at_t1 + E_x2_at_t2)
      
      if(log(runif(1))< -delta_energy){
        # Swap
        temp = current_x[idx1]
        current_x[idx1] = current_x[idx2]
        current_x[idx2] = temp
        swap_success = swap_success + 1
      }
      tot_iter = tot_iter + 1
    }
    chain_states[i, ] = current_x
  }
  cat("Swap Acceptance Rate:", round(swap_success/tot_iter, 3), "\n")
  return(list(chains = chain_states, t_ladder = t_ladder))
}

# Run and Visualize
res2 = run_PT_RW(iter = 5000, n_chains = 5)

# Visualize Chain 1 (The Target Chain)
par(mfrow=c(1,2))
plot(res2$chains[,1], type="l", col="blue", lwd=1,
     main="Chain 1 (t=0.01) - MYPT with RW", ylab="x", ylim = c(-3,3))
abline(h=c(-2, -1, 1, 2), col="red", lty=2)

# Visualize Histogram
hist(res2$chains[1:5000, 1], breaks=50, prob=TRUE, col="lightblue",
     main="Final Distribution (Parallel MY)", xlab="x", ylim = c(0,1.5),
     ylab = "Density")

# Overlay the True Target Density (Red Line)
curve(exp(-potential(x)), from = -4, to = 4, col = "red", lwd = 2, add = TRUE)

# Add legend
legend("topright", legend=c("PT Samples", "True Target"), 
       col=c("lightblue", "red"), lwd=c(10, 2))



# 3. Using MYPT with MALA
run_PT_MALA = function(iter, n_chains) {
  # T=0 is target. T=High is flat.
  t_ladder = c(0, exp(seq(log(0.01), log(5), length.out = n_chains-1)))
  
  # Step size
  tau = 0.05
  
  chain_states = matrix(0, nrow = iter, ncol = n_chains)
  current_x = rep(0, n_chains) # Start all at 0
  
  swap_success = 0
  tot_iter = 0
  
  for(i in 1:iter){
    for(k in 1:n_chains){
      t_val = t_ladder[k]
      prox_pt = moreau_func(current_x[k], t_val)$prox
      grad_old = ifelse (t_val != 0,  (1 / t_val) * (current_x[k] - prox_pt), 0)
      
      # Langevin Proposal
      mean_old = current_x[k] - 0.5 * tau * grad_old
      x_new = mean_old + sqrt(tau) * rnorm(1)
      
      # Calculate Energies
      U_old = moreau_func(current_x[k], t_val)$val
      U_new = moreau_func(x_new, t_val)$val
      
      # We need Gradient at x_new to calculate q factors
      prox_new = moreau_func(x_new, t_val)$prox
      grad_new = ifelse (t_val != 0, (1 / t_val) * (x_new - prox_new), 0)
      mean_new = x_new - 0.5 * tau * grad_new
      
      # Prob of proposing x_new given x_old
      log_q_fwd = dnorm(x_new, mean = mean_old, sd = sqrt(tau), log = TRUE)
      # Prob of proposing x_old given x_new (The Reverse Step)
      log_q_rev = dnorm(current_x[k], mean = mean_new, sd = sqrt(tau), log = TRUE)
      
      # MALA Acceptance Ratio
      log_ratio = (U_old - U_new) + (log_q_rev - log_q_fwd)
     
      if(log(runif(1)) < log_ratio){
        current_x[k] = x_new
      }
    }
    
    # Performing the swaps
    for (k in 1:(n_chains - 1)) {
      idx1 = k
      idx2 = k + 1
      
      x1 = current_x[idx1]
      x2 = current_x[idx2]
      t1 = t_ladder[idx1]
      t2 = t_ladder[idx2]
      
      # Ratio = exp(-[(M^t1(x2) + M^t2(x1)) - (M^t1(x1) + M^t2(x2))])
      
      # Calculate Energy of x2 if it moves to chain 1 (and vice versa)
      E_x2_at_t1 = moreau_func(x2, t1)$val
      E_x1_at_t2 = moreau_func(x1, t2)$val
      
      # Calculate Current Energies
      E_x1_at_t1 = moreau_func(x1, t1)$val
      E_x2_at_t2 = moreau_func(x2, t2)$val
      
      # Proposed Energy - Current Energy
      delta_energy = (E_x2_at_t1 + E_x1_at_t2) - (E_x1_at_t1 + E_x2_at_t2)
      
      if(log(runif(1))< -delta_energy){
        # Swap
        temp = current_x[idx1]
        current_x[idx1] = current_x[idx2]
        current_x[idx2] = temp
        swap_success = swap_success + 1
      }
      tot_iter = tot_iter + 1
    }
    chain_states[i, ] = current_x
  }
  cat("Swap Acceptance Rate:", round(swap_success/tot_iter, 3), "\n")
  return(list(chains = chain_states, t_ladder = t_ladder))
}

# Run and Visualize 
res3 = run_PT_MALA(iter = 5000, n_chains = 5)

# Visualize Chain 1 (Target T=1)
par(mfrow=c(1,2))
plot(res3$chains[,1], type="l", col="blue", lwd=1,
     main="Chain 1 (Target T=1) - MALA", ylab="x", ylim = c(-3,3))
abline(h=c(-2, -1, 1, 2), col="red", lty=2)

# Visualize Histogram
hist(res3$chains[1:5000, 1], breaks=50, prob=TRUE, col="lightblue",
     main="Final Distribution (PT with MALA)", xlab="x", ylim = c(0,1.5))
x_grid = seq(-4, 4, length=200)
lines(x_grid, exp(-potential(x_grid)), col="red", lwd=2)

# Add legend
legend("topright", legend=c("PT Samples", "True Target"), 
       col=c("lightblue", "red"), lwd=c(10, 2))


