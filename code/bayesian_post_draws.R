n <- 10
p <- 1
XtX <- 1
Xty <- 8
yty <- 65
lambda <- 2.4

x_range <- c(-1.5, 9)
y_range <- c(-4,4.5)
# Define the Log-Posterior in terms of (beta, ln(sigma^2))
# nu = ln(sigma^2)
log_posterior <- function(beta, nu) 
{
  # Quadratic term: (y - X*beta)^T (y - X*beta)
  Q <- yty - 2 * beta * Xty + (beta^2)* XtX
  # Log density (need to use Jacobian for change of variable)
  log_dens <- -((n-1)/2) * nu - (Q / (2 * exp(nu))) - lambda * abs(beta)
  return(log_dens)
}

# define a grid for calculations
beta_seq <- seq(-1.5, 9, length.out = 200)
nu_seq <- seq(-4, 4.5, length.out = 200)

# Bivariate Random Walk Metropolis-Hastings Algorithm with Parallel Tempering
run_rwmh <- function(iter, n_chains, init_beta, init_nu, step_beta, step_nu) 
{
  # T=1 is target. T=High is flat.
  temps = exp(seq(log(1), log(50), length.out = n_chains))
  betas = 1/temps 
  
  chain_beta = matrix(0, nrow = iter, ncol = n_chains)
  chain_nu = matrix(0, nrow = iter, ncol = n_chains)
  
  #Initialize
  curr_beta = rep(init_beta, n_chains)
  curr_nu = rep(init_nu, n_chains)
  
  tot_iter = 0
  swap_success = 0
  
  for(i in 1:iter){
    for(k in 1:n_chains){
      # MH Proposal
      prop_beta <- curr_beta[k] + rnorm(1, mean = 0, sd = step_beta)
      prop_nu <- curr_nu[k] + rnorm(1, mean = 0, sd = step_nu)
      
      # Evaluate Log Posteriors
      log_curr <- log_posterior(curr_beta[k], curr_nu[k])
      log_prop <- log_posterior(prop_beta, prop_nu)
      
      # Calculate Log Acceptance Ratio
      log_ratio <- betas[k]* (log_prop - log_curr)
      
      # Accept/Reject Step
      if(log(runif(1)) < log_ratio){
        curr_beta[k] <- prop_beta
        curr_nu[k] <- prop_nu
      }
    }
    # swapping chain k and k+1
    for (k in 1:(n_chains - 1)){
      chain_i = k
      chain_j = k + 1
      
      x_i = curr_beta[chain_i]
      x_j = curr_beta[chain_j]
      
      y_i = curr_nu[chain_i]
      y_j = curr_nu[chain_j]
      
      # The Swap Acceptance Formula = min(1, exp((Beta_i-Beta_j)*(U(x_i,y_i)-U(x_j,y_j))))
      d_beta = betas[chain_i]-betas[chain_j]
      d_energy = -log_posterior(x_i, y_i)+log_posterior(x_j, y_j)
      
      log_swap_ratio = d_beta * d_energy
      
      if(log(runif(1)) < log_swap_ratio){
        # Perform Swap
        temp_x = curr_beta[chain_i]
        curr_beta[chain_i] = curr_beta[chain_j]
        curr_beta[chain_j] = temp_x
        
        temp_y = curr_nu[chain_i]
        curr_nu[chain_i] = curr_nu[chain_j]
        curr_nu[chain_j] = temp_y
        
        swap_success = swap_success + 1
      }
      tot_iter = tot_iter + 1
    }
    # Store states
    chain_beta[i, ] = curr_beta
    chain_nu[i, ] = curr_nu
  }
  cat("Swap Acceptance Rate:", round(swap_success/tot_iter, 3), "\n")
  return(list(chain_beta = chain_beta, chain_nu = chain_nu, temps = temps))
}


# Intial values near the least sq estimates
res1 <- run_rwmh(iter = 20000,
                 n_chains = 6,
                    init_beta = 4, 
                    init_nu = 0, 
                    step_beta = 0.5,    
                    step_nu = 0.25)


# Plotting the Results
par(mfrow=c(1,2))

# Trace Plot for Beta
plot(res1$chain_beta[,1], type="l", col="blue", 
     main="Trace Plot of Beta", xlab="Iteration", ylab="Beta")
# Trace Plot for ln(sigma^2)
plot(res1$chain_nu[,1], type="l", col="blue", 
     main="Trace Plot of ln(Sigma2)", xlab="Iteration", ylab="Sigma2")

par(mfrow = c(1,1))
# Scatter plot between ln(sigma^2) and beta showing the sampled distribution
plot(res1$chain_beta[,1], res1$chain_nu[,1], 
     pch = 20, col = rgb(0, 0, 0, 0.1), 
     main = expression(paste("Sampled Joint Posterior of ", beta, " and ln(", sigma^2, ") using Random Walk")),
     xlab = expression(beta), ylab = expression(ln(sigma^2)),
     xlim = x_range, ylim = y_range)


#######################################################################
# Visualize Smoothed Moreau Contours for the t_ladder

plot_standard_pt_ladder = function(T_ladder, x_bounds, y_bounds, plot_grid_size)
{
  x_plot = seq(x_bounds[1], x_bounds[2], length.out = plot_grid_size)
  y_plot = seq(y_bounds[1], y_bounds[2], length.out = plot_grid_size)
  
  rows = floor(sqrt(length(T_ladder)))
  cols = ceiling(length(T_ladder) / rows)
  par(mfrow = c(rows, cols), mar = c(4, 4, 2, 1))
  
  for(Temp in T_ladder){
    cat("Calculating tempered contours for T =", round(Temp, 3), "...\n")
    
    log_density_matrix = matrix(0, nrow = plot_grid_size, ncol = plot_grid_size)
    
    for(i in 1:plot_grid_size) {
      for(j in 1:plot_grid_size) {
        # Standard PT: Divide the log-posterior by the Temperature
        log_density_matrix[i, j] = (1 / Temp) * log_posterior(x_plot[i], y_plot[j])
      }
    }
    
    # Stabilization Trick
    max_log_dens <- max(log_density_matrix)
    smooth_density <- exp(log_density_matrix - max_log_dens)
    
    # Set title based on whether it is the target or a tempered chain
    title_str = ifelse(Temp == 1, 
                       "Target Density (T = 1)", 
                       paste("Tempered Density (T =", round(Temp, 2), ")"))
    
    # Plot the tempered landscape
    contour(x_plot, y_plot, smooth_density, 
            levels = seq(0.01, 0.99, length.out = 25),
            col = "darkred", drawlabels = FALSE,
            main = title_str, 
            xlab = expression(beta), ylab = expression(ln(sigma^2)))
  }
  
  par(mfrow = c(1, 1))
  cat("Plotting complete!\n")
}

plot_standard_pt_ladder(T_ladder = res1$temps, 
                   x_bounds = c(min(beta_seq), max(beta_seq)), 
                   y_bounds = c(min(nu_seq), max(nu_seq)), 
                   plot_grid_size = 100)

##########################################


# Bivariate Random Walk Metropolis-Hastings Algorithm using Moreau Envelopes approximation

# Potential is the -1*log_posterior
potential <- function(beta, nu){
  return(-log_posterior(beta,nu))
}

# define a grid for calculations
beta_seq <- seq(-1.5, 9, length.out = 200)
nu_seq <- seq(-4, 4.5, length.out = 200)

# Pre-calculate the Potential on this grid
potential_on_grid <- outer(beta_seq, nu_seq, Vectorize(potential))

# Moreau Envelope calculation
moreau_func_2d <- function(beta, nu, t) 
{
  if(t==0){
    return(potential(beta,nu))
  }
  penalty <- outer(beta_seq, nu_seq, function(b, n_val) {
    (1/(2*t)) * ((beta - b)^2 + (nu- n_val)^2)
  })
  objectives <- potential_on_grid + penalty
  min_value <- min(objectives)                
  return(min_value)
}


# Standard RW MH on the Moreau Envelope
run_PTMY_RW <- function(iter, n_chains, init_beta, init_nu, step_beta, step_nu)
{
  # T=0 is target. T=High is flat.
  t_ladder = c(0, exp(seq(log(0.01), log(100), length.out = n_chains-1)))
  
  # Step size
  tau = 0.05
  
  chain_beta = matrix(0, nrow = iter, ncol = n_chains)
  chain_nu = matrix(0, nrow = iter, ncol = n_chains)
  
  #Initialize
  curr_beta = rep(init_beta, n_chains)
  curr_nu = rep(init_nu, n_chains)
  
  tot_iter = 0
  swap_success = 0
  
  for (i in 1:iter) {
    for(k in 1:n_chains){
      t_val = t_ladder[k]
      
      # MH Proposal
      prop_beta <- curr_beta[k] + rnorm(1, mean = 0, sd = step_beta)
      prop_nu <- curr_nu[k] + rnorm(1, mean = 0, sd = step_nu)
      
      # Calculate Moreau Energy 
      U_old <- moreau_func_2d(curr_beta[k], curr_nu[k], t_val)
      U_new <- moreau_func_2d(prop_beta, prop_nu, t_val)
      
      # Acceptance Ratio
      log_ratio <- (U_old - U_new)
      
      # Accept/Reject Step
      if(log(runif(1)) < log_ratio){
        curr_beta[k] <- prop_beta
        curr_nu[k] <- prop_nu
      }
    }
    
    # Performing the swaps
    for(k in 1:(n_chains - 1)){
      idx1 = k
      idx2 = k + 1
      
      x1 = curr_beta[idx1]
      x2 = curr_beta[idx2]
      
      y1 = curr_nu[idx1]
      y2 = curr_nu[idx2]
      
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
        temp_x = curr_beta[idx1]
        curr_beta[idx1] = curr_beta[idx2]
        curr_beta[idx2] = temp_x
        
        temp_y = curr_nu[idx1]
        curr_nu[idx1] = curr_nu[idx2]
        curr_nu[idx2] = temp_y
        
        swap_success = swap_success + 1
      }
      tot_iter = tot_iter + 1
    }
    chain_beta[i, ] = curr_beta
    chain_nu[i, ] = curr_nu
  }
  cat("Swap Acceptance Rate:", round(swap_success/tot_iter, 3), "\n")
  return(list(chain_beta = chain_beta, chain_nu = chain_nu, t_ladder = t_ladder))
}



# Intial values near the least sq estimates
res2 <- run_PTMY_RW(iter = 20000,
                 n_chains = 6,
                 init_beta = 4, 
                 init_nu = 0, 
                 step_beta = 0.5,    
                 step_nu = 0.25)


# Plotting the Results
par(mfrow=c(1,2))

# Trace Plot for Beta
plot(res2$chain_beta[,1], type="l", col="blue", 
     main="Trace Plot of Beta", xlab="Iteration", ylab="Beta")
# Trace Plot for ln(sigma^2)
plot(res2$chain_nu[,1], type="l", col="blue", 
     main="Trace Plot of ln(Sigma2)", xlab="Iteration", ylab="Sigma2")

par(mfrow = c(1,1))
# Scatter plot between ln(sigma^2) and beta showing the sampled distribution
plot(res2$chain_beta[,1], res2$chain_nu[,1], 
     pch = 20, col = rgb(0, 0, 0, 0.1), 
     main = expression(paste("Sampled Joint Posterior of ", beta, " and ln(", sigma^2, ") using Moreau RW")),
     xlab = expression(beta), ylab = expression(ln(sigma^2)),
     xlim = x_range, ylim = y_range)


#############################################################


# Visualize Smoothed Moreau Contours for the t_ladder

plot_moreau_ladder = function(t_ladder, x_bounds, y_bounds, plot_grid_size)
{
  x_plot = seq(x_bounds[1], x_bounds[2], length.out = plot_grid_size)
  y_plot = seq(y_bounds[1], y_bounds[2], length.out = plot_grid_size)
  
  rows = floor(sqrt(length(t_ladder)))
  cols = ceiling(length(t_ladder) / rows)
  par(mfrow = c(rows, cols), mar = c(4, 4, 2, 1))
  
  for(t in t_ladder){
    if(t == 0){
      # Calculate log posterior over the grid
      z_log <- outer(beta_seq, nu_seq, Vectorize(log_posterior))
      
      # To avoid underflow/overflow when exponentiating, subtract the maximum log-density
      z_log_max <- max(z_log)
      z <- exp(z_log - z_log_max)
      
      # t=0 is the original target density
      contour(beta_seq, nu_seq, z, 
              levels = seq(0.01, 0.99, length.out = 25),
              col = "darkred", drawlabels = FALSE,
              main = "Target Density (t = 0)", xlab = expression(beta), ylab = expression(ln(sigma^2)))
    } 
    else{
      cat("Calculating smoothed contours for t =", round(t, 3), "...\n")
      
      log_density_matrix = matrix(0, nrow = plot_grid_size, ncol = plot_grid_size)
      
      for(i in 1:plot_grid_size) {
        for(j in 1:plot_grid_size) {
          smoothed_U = moreau_func_2d(x_plot[i], y_plot[j], t)
          log_density_matrix[i, j] = -smoothed_U
        }
      }
      
      max_log_dens <- max(log_density_matrix)
      smooth_density <- exp(log_density_matrix - max_log_dens)
      
      # Plot the smoothed landscape
      contour(x_plot, y_plot, smooth_density, 
              levels = seq(0.01, 0.99, length.out = 25),
              col = "darkred", drawlabels = FALSE,
              main = paste("Smoothed Density (t =", round(t, 2), ")"), 
              xlab = expression(beta), ylab = expression(ln(sigma^2)))
    }
  }
  
  par(mfrow = c(1, 1))
  cat("Plotting complete!\n")
}

plot_moreau_ladder(t_ladder = res2$t_ladder, 
                   x_bounds = c(min(beta_seq), max(beta_seq)), 
                   y_bounds = c(min(nu_seq), max(nu_seq)), 
                   plot_grid_size = 100)

##########################################

# U = F + G , G = lambda|beta|
# The Analytical Moreau Envelope for the lambda|beta| - L1 Penalty (Huber Function)
moreau_penalty_analytical <- function(beta, lambda, t)
{
  if (t == 0) {
    return(lambda * abs(beta))
  }
  threshold <- t * lambda
  if (abs(beta) <= threshold) {
    return((beta^2) / (2 * t))
  } else {
    return(lambda * abs(beta) - (t * lambda^2) / 2)
  }
}

# The Total Analytical Smoothed Potential (F + M_g)
analytical_moreau_2d <- function(beta, nu, t)
{
  Q <- yty - 2 * beta * Xty + (beta^2) * XtX
  f <- ((n-1)/2) * nu + (Q / (2 * exp(nu)))
  moreau_penalty <- moreau_penalty_analytical(beta, lambda, t)
  return(f + moreau_penalty)
}

run_PTMY_RW_analytical <- function(iter, n_chains, init_beta, init_nu, step_beta, step_nu)
{
  # T=0 is target. T=High is flat.
  t_ladder = c(0, exp(seq(log(0.01), log(5), length.out = n_chains-1)))
  
  # Step size
  tau = 0.05
  
  chain_beta = matrix(0, nrow = iter, ncol = n_chains)
  chain_nu = matrix(0, nrow = iter, ncol = n_chains)
  
  #Initialize
  curr_beta = rep(init_beta, n_chains)
  curr_nu = rep(init_nu, n_chains)
  
  tot_iter = 0
  swap_success = 0
  
  for (i in 1:iter) {
    for(k in 1:n_chains){
      t_val = t_ladder[k]
      
      # MH Proposal
      prop_beta <- curr_beta[k] + rnorm(1, mean = 0, sd = step_beta)
      prop_nu <- curr_nu[k] + rnorm(1, mean = 0, sd = step_nu)
      
      # Calculate Moreau Energy 
      U_old <- analytical_moreau_2d(curr_beta[k], curr_nu[k], t_val)
      U_new <- analytical_moreau_2d(prop_beta, prop_nu, t_val)
      
      # Acceptance Ratio
      log_ratio <- (U_old - U_new)
      
      # Accept/Reject Step
      if(log(runif(1)) < log_ratio){
        curr_beta[k] <- prop_beta
        curr_nu[k] <- prop_nu
      }
    }
    
    # Performing the swaps
    for(k in 1:(n_chains - 1)){
      idx1 = k
      idx2 = k + 1
      
      x1 = curr_beta[idx1]
      x2 = curr_beta[idx2]
      
      y1 = curr_nu[idx1]
      y2 = curr_nu[idx2]
      
      t1 = t_ladder[idx1]
      t2 = t_ladder[idx2]
      
      # Ratio = exp(-[(M^t1(x2,y2) + M^t2(x1,y1)) - (M^t1(x1,y1) + M^t2(x2,y2))])
      
      # Calculate Energy of (x2,y2) if it moves to chain 1 (and vice versa)
      E2_at_t1 = analytical_moreau_2d(x2, y2, t1)
      E1_at_t2 = analytical_moreau_2d(x1, y1, t2)
      
      # Calculate Current Energies
      E1_at_t1 = analytical_moreau_2d(x1, y1, t1)
      E2_at_t2 = analytical_moreau_2d(x2, y2, t2)
      
      # Proposed Energy - Current Energy
      delta_energy = (E2_at_t1 + E1_at_t2) - (E1_at_t1 + E2_at_t2)
      
      if(log(runif(1))< -delta_energy){
        # Swap
        temp_x = curr_beta[idx1]
        curr_beta[idx1] = curr_beta[idx2]
        curr_beta[idx2] = temp_x
        
        temp_y = curr_nu[idx1]
        curr_nu[idx1] = curr_nu[idx2]
        curr_nu[idx2] = temp_y
        
        swap_success = swap_success + 1
      }
      tot_iter = tot_iter + 1
    }
    chain_beta[i, ] = curr_beta
    chain_nu[i, ] = curr_nu
  }
  cat("Swap Acceptance Rate:", round(swap_success/tot_iter, 3), "\n")
  return(list(chain_beta = chain_beta, chain_nu = chain_nu, t_ladder = t_ladder))
}



# Intial values near the least sq estimates
res3 <- run_PTMY_RW_analytical(iter = 20000,
                    n_chains = 6,
                    init_beta = 4, 
                    init_nu = 0, 
                    step_beta = 0.5,    
                    step_nu = 0.25)


# Plotting the Results
par(mfrow=c(1,2))

# Trace Plot for Beta
plot(res3$chain_beta[,1], type="l", col="blue", 
     main="Trace Plot of Beta", xlab="Iteration", ylab="Beta")
# Trace Plot for ln(sigma^2)
plot(res3$chain_nu[,1], type="l", col="blue", 
     main="Trace Plot of ln(Sigma2)", xlab="Iteration", ylab="Sigma2")

par(mfrow = c(1,1))
# Scatter plot between ln(sigma^2) and beta showing the sampled distribution
plot(res3$chain_beta[,1], res3$chain_nu[,1], 
     pch = 20, col = rgb(0, 0, 0, 0.1), 
     main = expression(paste("Sampled Joint Posterior of ", beta, " and ln(", sigma^2, ") using analytical Moreau RW")),
     xlab = expression(beta), ylab = expression(ln(sigma^2)),
     xlim = x_range, ylim = y_range)


###################################


# Visualize Smoothed Moreau Contours for the t_ladder

plot_moreau_ladder = function(t_ladder, x_bounds, y_bounds, plot_grid_size)
{
  x_plot = seq(x_bounds[1], x_bounds[2], length.out = plot_grid_size)
  y_plot = seq(y_bounds[1], y_bounds[2], length.out = plot_grid_size)
  
  rows = floor(sqrt(length(t_ladder)))
  cols = ceiling(length(t_ladder) / rows)
  par(mfrow = c(rows, cols), mar = c(4, 4, 2, 1))
  
  for(t in t_ladder){
    if(t == 0){
      # Calculate log posterior over the grid
      z_log <- outer(beta_seq, nu_seq, Vectorize(log_posterior))
      
      # To avoid underflow/overflow when exponentiating, subtract the maximum log-density
      z_log_max <- max(z_log)
      z <- exp(z_log - z_log_max)
      
      # t=0 is the original target density
      contour(beta_seq, nu_seq, z, 
              levels = seq(0.01, 0.99, length.out = 25),
              col = "darkred", drawlabels = FALSE,
              main = "Target Density (t = 0)", xlab = expression(beta), ylab = expression(ln(sigma^2)))
    } 
    else{
      cat("Calculating smoothed contours for t =", round(t, 3), "...\n")
      
      log_density_matrix = matrix(0, nrow = plot_grid_size, ncol = plot_grid_size)
  
      for(i in 1:plot_grid_size) {
        for(j in 1:plot_grid_size) {
          smoothed_U = analytical_moreau_2d(x_plot[i], y_plot[j], t)
          log_density_matrix[i, j] = -smoothed_U
        }
      }
      
      max_log_dens <- max(log_density_matrix)
      smooth_density <- exp(log_density_matrix - max_log_dens)
      
      # Plot the smoothed landscape
      contour(x_plot, y_plot, smooth_density, 
              levels = seq(0.01, 0.99, length.out = 25),
              col = "darkred", drawlabels = FALSE,
              main = paste("Smoothed Density (t =", round(t, 2), ")"), 
              xlab = expression(beta), ylab = expression(ln(sigma^2)))
    }
  }
  
  par(mfrow = c(1, 1))
  cat("Plotting complete!\n")
}

plot_moreau_ladder(t_ladder = res3$t_ladder, 
                   x_bounds = c(min(beta_seq), max(beta_seq)), 
                   y_bounds = c(min(nu_seq), max(nu_seq)), 
                   plot_grid_size = 100)

############################

# expected squared jump distance
esjd = function(chain_x, chain_y)
{
  chain = cbind(chain_x, chain_y)
  diffs = diff(chain)
  sq_dists = rowSums(diffs^2)
  return(mean(sq_dists))
}

n_reps = 10

esjd_rw = numeric(n_reps)
esjd_moreau = numeric(n_reps)
esjd_moreau_anal = numeric(n_reps)

for(i in 1:n_reps) {
  res1 <- run_rwmh(iter = 5000,
                   n_chains = 6,
                   init_beta = 4, 
                   init_nu = 0, 
                   step_beta = 0.5,    
                   step_nu = 0.25)
  res2 <- run_PTMY_RW(iter = 5000,
                      n_chains = 6,
                      init_beta = 4, 
                      init_nu = 0, 
                      step_beta = 0.5,    
                      step_nu = 0.25)
  res3 <- run_PTMY_RW_analytical(iter = 5000,
                                 n_chains = 6,
                                 init_beta = 4, 
                                 init_nu = 0, 
                                 step_beta = 0.5,    
                                 step_nu = 0.25)
  
  esjd_rw[i] = esjd(res1$chain_beta[,1], res1$chain_nu[,1])
  esjd_moreau[i] = esjd(res2$chain_beta[,1], res2$chain_nu[,1])
  esjd_moreau_anal[i] = esjd(res3$chain_beta[,1], res3$chain_nu[,1])
}

avg_esjd_pt = mean(esjd_rw)
avg_esjd_mypt = mean(esjd_moreau)
avg_esjd_mypt_anal = mean(esjd_moreau_anal)
cat("Average ESJD (Standard PT):", round(avg_esjd_pt, 3), "\n")
cat("Average ESJD (MYPT):", round(avg_esjd_mypt, 3), "\n")
cat("Average ESJD (Analytical MYPT):", round(avg_esjd_mypt_anal, 3), "\n")

#######################################
