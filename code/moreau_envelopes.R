# Moreau envelopes of a Gaussian mixture
 
# Target Potential U(x) is a Mixture of Gaussians
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

# Moreau Envelope Calculation
# we solve the optimization problem: min_y { U(y) + (1/2t) * ||x - y||^2 }
calc_moreau_val = function(x, t){
  if(t==0){
    return(potential(x))
  }
  penalty = (1/(2*t))*(x-grid)^2
  all_objectives = potential_on_grid + penalty
  return(min(all_objectives)) # Return the minimized value
}

density_val = function(x, t){
  return(exp(-1*calc_moreau_val(x, t)))
}

# Plotting
# Plot for different values of t
t_values = c(0, 0.01, 0.1, 0.5, 2, 10, 100) 
colors = c("black", "blue", "lightblue", "gold", "orange", "red", "maroon")

# Potential Energy plot
par(mfrow=c(1,2))
plot(NULL, xlim=c(-4, 4), ylim=c(0, 11), 
     xlab="x", ylab="Potential Energy", main="Moreau Envelopes")

for (i in 1:length(t_values)){
  t = t_values[i]
  if(t==0){
    y_seq = potential_on_grid
  } 
  else{
    y_seq = sapply(grid, function(x) calc_moreau_val(x, t))
  }
  lines(grid, y_seq, col=colors[i], lwd=2)
}
legend("top", legend=paste("t =", t_values), col=colors, lwd=2, bty="n", cex = 0.5)


# Density plot
plot(NULL, xlim=c(-4, 4), ylim=c(0, 2.1), 
     xlab="x", ylab="Density", main="Moreau Envelopes")

for (i in 1:length(t_values)){
  t = t_values[i]
  if(t==0){
    y_seq = exp(-potential_on_grid)
  }
  else{
    y_seq = sapply(grid, function(x) density_val(x, t))
  }
  lines(grid, y_seq, col=colors[i], lwd=2)
}
legend("top", legend=paste("t =", t_values), col=colors, lwd=2, bty="n", cex = 0.5)


# Calculate Proximal Operator
# This returns the optimal 'y' that minimizes: G(y) + (1/2t)||x-y||^2
get_prox = function(x_current, t){
  if(t==0) return(x_current)
  penalty = (1/(2*t))*(x_current - grid)^2
  objectives = potential_on_grid + penalty
  min_idx = which.min(objectives)
  return(grid[min_idx])
}

# DAZ Algorithm Implementation
run_daz = function(initial_state, n_levels = 50, K = 20) {
  t_schedule =  exp(seq(log(.01), log(0.0001), length.out = n_levels))
  tau_schedule = 0.5 * t_schedule
  
  # Initialize
  current_x = initial_state
  trajectory = numeric(n_levels * K)
  step_counter = 1

  for (n in 1:n_levels) {
    t_val = t_schedule[n]
    tau_val = tau_schedule[n]
    
    for (k in 1:K) {
      prox_pt = get_prox(current_x, t_val)
      grad_moreau = (1 / t_val) * (current_x - prox_pt)
      noise = rnorm(1)
      current_x = current_x - (tau_val * grad_moreau) + sqrt(2 * tau_val) * noise
      trajectory[step_counter] = current_x
      step_counter = step_counter + 1
    }
  }
  return(list(traj = trajectory, final = current_x))
}

# Run & Visualize
#set.seed(123)

# Start at x=0 
res = run_daz(initial_state = 0,
               n_levels = 50,
               K = 20)

# Plot the path
# X-axis: Time steps
# Y-axis: Position of the walker
par(mfrow = c(1,2))
plot(res$traj, type="l", col="darkblue", lwd=1,
     main="DAZ Sampling Trajectory", xlab="Time Steps", ylab="Position (x)",
     ylim=c(-3, 3))

# Add lines for the true modes
abline(h=c(-2, -1, 1, 2), col="red", lty=2)
text(0, 2.2, "True Modes", col="red", pos=4)
text(0, res$traj[1], "Start (0)", pos=4)

plot(density(res$traj), main="Density of DAZ Samples", xlab="x", ylab="Density",
     col="darkgreen", lwd=2, ylim = c(0,1.5))
# Overlay true density
true_density = sapply(grid, potential)
lines(grid, exp(-true_density), col="red", lwd=2)
legend("topleft", legend=c("DAZ Samples", "True Density"),
       col=c("darkgreen", "red"), lwd=2, bty="n")



# n = 5000, k = 2000
# Start at x = 0 
res2 = run_daz(initial_state = 0,
              n_levels = 5000,
              K = 2000)

# Plot the path
# X-axis: Time steps
# Y-axis: Position of the walker
par(mfrow = c(1,2))
samples = res2$traj[seq(1, length(res$traj), by = 1000)]
plot(samples, type="l", col="darkblue", lwd=1,
     main="DAZ Sampling Trajectory", xlab="Time Steps", ylab="Position (x)",
     ylim=c(-3, 3))

# Add lines for the true modes
abline(h=c(-2, -1, 1, 2), col="red", lty=2)
text(0, 2.2, "True Modes", col="red", pos=4)
text(0, res$traj[1], "Start (0)", pos=4)

plot(density(samples), main="Density of DAZ Samples", xlab="x", ylab="Density",
     col="darkgreen", lwd=2, ylim = c(0,1.5))
# Overlay true density
true_density = sapply(grid, potential)
lines(grid, exp(-true_density), col="red", lwd=2)
legend("topleft", legend=c("DAZ Samples", "True Density"),
       col=c("darkgreen", "red"), lwd=2, bty="n")


