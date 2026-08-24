# Bivariate gaussian with 20 modes
set.seed(123) 
n_modes = 20

# Randomize Means (between 0 and 10)
mu_x = runif(n_modes, min = 0, max = 10)
mu_y = runif(n_modes, min = 0, max = 8)

# Randomize Covariances
sigma_x = runif(n_modes, min = 0.2, max = 0.5)
sigma_y = runif(n_modes, min = 0.2, max = 0.4)
rho = runif(n_modes, min = -0.4, max = 0.4) 

# Randomize Weights (Normalize so they sum to 1)
raw_weights = runif(n_modes, min = 0.4, max = 0.7)
weights = raw_weights / sum(raw_weights)

# bivariate normal density
biv_norm = function(x, y, mx, my, sx, sy, r) 
{
  z = ((x - mx)^2 / sx^2) + ((y - my)^2 / sy^2) - (2 * r * (x - mx) * (y - my) / (sx * sy))
  coef = 1 / (2 * pi * sx * sy * sqrt(1 - r^2)) # normalizing constant
  return(coef * exp(-z / (2 * (1 - r^2))))
}

# The Mixture Density Function
density_func = function(x, y) 
{
  total_density = 0
  # Sum the weighted density of all 20 components
  for (i in 1:n_modes) {
    comp_density = biv_norm(x, y, mu_x[i], mu_y[i], sigma_x[i], sigma_y[i], rho[i])
    total_density = total_density + weights[i] * comp_density
  }
  return(total_density)
}

# Vectorize
density_func_vec = Vectorize(density_func)

# Plotting the density
x_seq = seq(-1, 11, length.out = 200)
y_seq = seq(-1, 9, length.out = 200)

z_density = outer(x_seq, y_seq, density_func_vec)

contour(x_seq, y_seq, z_density, 
        nlevels = 20, col = "darkblue", drawlabels = FALSE,
        main = "Multi-Modal Mixture Gaussian Density",
        xlab = "X", ylab = "Y")

potential = function(x, y) 
{
  dens_val = density_func(x, y)
  # Add 1e-300 to prevent log(0)
  return(-log(dens_val + 1e-300)) 
}

#######################

# Samples from distribution using Parallel Tempering with Random Walk
run_PT_RW = function(iter, n_chains, init_x, init_y, sigma_x, sigma_y) 
{
  # T=1 is target. T=High is flat.
  temps = exp(seq(log(1), log(50), length.out = n_chains))
  betas = 1/temps 
  
  chain_states_x = matrix(0, nrow = iter, ncol = n_chains)
  chain_states_y = matrix(0, nrow = iter, ncol = n_chains)
  
  #Initialize
  current_x = rep(init_x, n_chains)
  current_y = rep(init_y, n_chains)
  
  tot_iter = 0
  swap_success = 0
  acceptance = numeric(length = n_chains)
  for(i in 1:iter){
    for(k in 1:n_chains){
      # Propose new value
      x_old = current_x[k]
      x_new = x_old + rnorm(1, mean=0, sd= sqrt(temps[k])* sigma_x) # Random walk
      
      y_old = current_y[k]
      y_new = y_old + rnorm(1, mean=0, sd=sqrt(temps[k])* sigma_y) # Random walk
      
      # Calculate Energy Difference
      U_old = potential(x_old, y_old)
      U_new = potential(x_new, y_new)
      
      # Acceptance Ratio = exp(-beta *(U_new - U_old))
      log_ratio = -betas[k] *(U_new - U_old)
      
      if(log(runif(1)) < log_ratio){
        current_x[k] = x_new
        current_y[k] = y_new
        acceptance[k] = acceptance[k] + 1
      }
    }
    
    # swapping chain k and k+1
    for (k in 1:(n_chains - 1)){
      chain_i = k
      chain_j = k + 1
      
      x_i = current_x[chain_i]
      x_j = current_x[chain_j]
      
      y_i = current_y[chain_i]
      y_j = current_y[chain_j]
      
      # The Swap Acceptance Formula = min(1, exp((Beta_i-Beta_j)*(U(x_i,y_i)-U(x_j,y_j))))
      d_beta = betas[chain_i]-betas[chain_j]
      d_energy = potential(x_i, y_i)-potential(x_j, y_j)
      
      log_swap_ratio = d_beta * d_energy
      
      if(log(runif(1)) < log_swap_ratio){
        # Perform Swap
        temp_x = current_x[chain_i]
        current_x[chain_i] = current_x[chain_j]
        current_x[chain_j] = temp_x
        
        temp_y = current_y[chain_i]
        current_y[chain_i] = current_y[chain_j]
        current_y[chain_j] = temp_y
        
        swap_success = swap_success + 1
      }
      tot_iter = tot_iter + 1
    }
    # Store states
    chain_states_x[i, ] = current_x
    chain_states_y[i, ] = current_y
  }
  
  cat("Swap Acceptance Rate:", round(swap_success/tot_iter, 3), "\n")
  cat("RW Acceptance Rates per Chain:", round(acceptance/iter, 3), "\n")
  return(list(chain_x = chain_states_x, chain_y = chain_states_y, temps = temps))
}


# Run and Visualize
res1 = run_PT_RW(iter = 5000, n_chains = 5, 
                       init_x = 5, init_y = 6, 
                       sigma_x = 0.4, sigma_y = 0.4)

# Draw the contour plot and the samples on the same plot
contour(x_seq, y_seq, z_density, 
        nlevels = 20, col = "darkblue", drawlabels = FALSE,
        main = "Samples using Random Walk",
        xlab = "X", ylab = "Y")
points(res1$chain_x[,1], res1$chain_y[,1], pch = 20, 
       col = rgb(0.8, 0, 0.4, 0.3), cex = 0.7)


#####################

# Samples from distribution using Parallel Tempering with Moreau Envelopes Random Walk

# Pre-calculate the Potential on this grid
potential_on_grid = outer(x_seq, y_seq, Vectorize(potential))

# Moreau Envelope calculation
moreau_func_2d = function(x_val, y_val, t) 
{
  if(t == 0){
    return(potential(x_val, y_val))
  }
  penalty = outer(x_seq, y_seq, function(x, y) {
    (1/(2*t)) * ((x_val - x)^2 + (y_val - y)^2)
  })
  objectives = potential_on_grid + penalty
  min_value = min(objectives)        
  return(min_value) 
}

run_MYPT_RW = function(iter, n_chains, init_x, init_y, sigma_x, sigma_y) 
{
  # T=0 is target. T=High is flat.
  t_ladder = c(0, exp(seq(log(0.01), log(100), length.out = n_chains-1)))
  
  # Step size
  tau = 0.05
  
  chain_states_x = matrix(0, nrow = iter, ncol = n_chains)
  chain_states_y = matrix(0, nrow = iter, ncol = n_chains)
  
  # Initialize
  current_x = rep(init_x, n_chains) 
  current_y = rep(init_y, n_chains)
  
  swap_success = 0
  tot_iter = 0
  acceptance = numeric(length = n_chains)
  for(i in 1:iter){
    for(k in 1:n_chains){
      t_val = t_ladder[k]
      
      # RW proposal (symmetric) 
      x_old = current_x[k]
      x_new = x_old + rnorm(1, mean=0, sd= sqrt(1 + 5*t_val)*sigma_x)
      y_old = current_y[k]
      y_new = y_old + rnorm(1, mean=0, sd=sqrt(1 + 5*t_val)*sigma_y) 
      
      # RW accept-reject step
      # Calculate Energy Difference
      U_old = moreau_func_2d(x_old, y_old, t_val)
      U_new = moreau_func_2d(x_new, y_new, t_val)
      
      # Acceptance Ratio
      log_ratio = (U_old - U_new)
      
      if(log(runif(1)) < log_ratio){
        current_x[k] = x_new
        current_y[k] = y_new
        acceptance[k] = acceptance[k] + 1
      }
    }
    
    # Performing the swaps
    for(k in 1:(n_chains - 1)){
      idx1 = k
      idx2 = k + 1
      
      x1 = current_x[idx1]
      x2 = current_x[idx2]
      
      y1 = current_y[idx1]
      y2 = current_y[idx2]
      
      t1 = t_ladder[idx1]
      t2 = t_ladder[idx2]
      
      # Ratio = exp(-[(M^t1(x2,y2) + M^t2(x1,y1)) - (M^t1(x1,y1) + M^t2(x2,y2))])
  
      # Calculate Energy of (x2,y2) if it moves to chain 1 (and vice versa)
      E2_at_t1 = moreau_func_2d(x2, y2, t1)
      E1_at_t2 = moreau_func_2d(x1, y1, t2)
      
      # Calculate Current Energies
      E1_at_t1 = moreau_func_2d(x1, y1, t1)
      E2_at_t2 = moreau_func_2d(x2, y2, t2)
      
      # Proposed Energy - Current Energy
      delta_energy = (E2_at_t1 + E1_at_t2) - (E1_at_t1 + E2_at_t2)
      
      if(log(runif(1))< -delta_energy){
        # Swap
        temp_x = current_x[idx1]
        current_x[idx1] = current_x[idx2]
        current_x[idx2] = temp_x
        
        temp_y = current_y[idx1]
        current_y[idx1] = current_y[idx2]
        current_y[idx2] = temp_y
        
        swap_success = swap_success + 1
      }
      tot_iter = tot_iter + 1
    }
    chain_states_x[i, ] = current_x
    chain_states_y[i, ] = current_y
  }
  cat("Swap Acceptance Rate:", round(swap_success/tot_iter, 3), "\n")
  cat("RW Acceptance Rates per Chain:", round(acceptance/iter, 3), "\n")
  return(list(chain_x = chain_states_x, chain_y = chain_states_y, t_ladder = t_ladder))
}

# Run and Visualize
res2 = run_MYPT_RW(iter = 5000, n_chains = 6, 
                 init_x = 5, init_y = 6, 
                 sigma_x = 0.4, sigma_y = 0.4)

# Draw the contour and the samples on the same plot
contour(x_seq, y_seq, z_density, 
        nlevels = 20, col = "darkblue", drawlabels = FALSE,
        main = "Samples using Moreau Envelopes RW",
        xlab = "X", ylab = "Y")
points(res2$chain_x[,1], res2$chain_y[,1], pch = 20, 
       col = rgb(0.8, 0, 0.4, 0.3), cex = 0.7)


#########################################################

# Visualize Smoothed Moreau Contours for the t_ladder
plot_moreau_ladder = function(t_ladder, x_bounds, y_bounds, plot_grid_size)
{
  # Use a coarser grid for visualization to keep computation fast
  x_plot = seq(x_bounds[1], x_bounds[2], length.out = plot_grid_size)
  y_plot = seq(y_bounds[1], y_bounds[2], length.out = plot_grid_size)
  
  # Set up a multi-panel plot based on the number of chains
  rows = floor(sqrt(length(t_ladder)))
  cols = ceiling(length(t_ladder) / rows)
  par(mfrow = c(rows, cols), mar = c(4, 4, 2, 1))
  
  for(t in t_ladder){
    if(t == 0){
      # t=0 is just the original target density
      contour(x_seq, y_seq, z_density, 
              nlevels = 20, col = "darkred", drawlabels = FALSE,
              main = "Target Density (t = 0)", xlab = "X", ylab = "Y")
    } 
    else{
      cat("Calculating smoothed contours for t =", round(t, 3), "...\n")
      
      # Initialize matrix for smoothed density
      smooth_density = matrix(0, nrow = plot_grid_size, ncol = plot_grid_size)
      
      # Calculate Moreau envelope for each point on the plotting grid
      for(i in 1:plot_grid_size) {
        for(j in 1:plot_grid_size) {
          # Get smoothed potential and transform to density: exp(-U)
          smoothed_U = moreau_func_2d(x_plot[i], y_plot[j], t)
          smooth_density[i, j] = exp(-smoothed_U)
        }
      }
      
      # Plot the smoothed landscape
      contour(x_plot, y_plot, smooth_density, 
              nlevels = 20, col = "darkred", drawlabels = FALSE,
              main = paste("Smoothed Density (t =", round(t, 2), ")"), 
              xlab = "X", ylab = "Y")
    }
  }
  
  # Reset plot layout
  par(mfrow = c(1, 1))
  cat("Plotting complete!\n")
}

# Run the visualization using the t_ladder from your results
# We pass the min/max of your original x_seq and y_seq
plot_moreau_ladder(t_ladder = res2$t_ladder, 
                   x_bounds = c(min(x_seq), max(x_seq)), 
                   y_bounds = c(min(y_seq), max(y_seq)), 
                   plot_grid_size = 100) 

##########################

# expected squared jump distance
esjd = function(chain_x, chain_y)
{
  chain = cbind(chain_x, chain_y)
  diffs = diff(chain)
  sq_dists = rowSums(diffs^2)
  return(mean(sq_dists))
}

n_reps = 100

esjd_rw = numeric(n_reps)
esjd_moreau = numeric(n_reps)

for(i in 1:n_reps) {
  res1 = run_PT_RW(iter = 5000, n_chains = 5, 
                   init_x = 5, init_y = 6, 
                   sigma_x = 0.4, sigma_y = 0.4)
  res2 = run_MYPT_RW(iter = 5000, n_chains = 6, 
                     init_x = 5, init_y = 6, 
                     sigma_x = 0.4, sigma_y = 0.4)

  esjd_rw[i] = esjd(res1$chain_x[,1], res1$chain_y[,1])
  esjd_moreau[i] = esjd(res2$chain_x[,1], res2$chain_y[,1])
}

avg_esjd_pt = mean(esjd_rw)
avg_esjd_mypt = mean(esjd_moreau)

cat("Average ESJD (Standard PT):", round(avg_esjd_pt, 3), "\n")
cat("Average ESJD (MYPT):", round(avg_esjd_mypt, 3), "\n")

#############################