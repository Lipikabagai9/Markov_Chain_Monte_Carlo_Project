# Define the grid ranges of beta and ln(sigma^2)
# nu = ln(sigma^2)
beta_seq <- seq(-1.5, 9, length.out = 300)
nu_seq <- seq(-4, 4.5, length.out = 300)

n <- 10
p <- 1
XtX <- 1
Xty <- 8
yty <- 65
lambda <- 2.4

# Function to calculate the log posterior for (beta, ln(sigma^2))
log_posterior <- function(beta, nu) {
  # Quadratic term: (y - X*beta)^T (y - X*beta)
  Q <- yty - 2 * beta * Xty + (beta^2)* XtX
  
  # Log density (need to use Jacobian for change of variable)
  log_dens <- -((n-1)/2) * nu - (Q / (2 * exp(nu))) - lambda * abs(beta)
  return(log_dens)
}

# Calculate log posterior over the grid
z_log <- outer(beta_seq, nu_seq, Vectorize(log_posterior))

# To avoid underflow/overflow when exponentiating, subtract the maximum log-density
z_log_max <- max(z_log)
z <- exp(z_log - z_log_max)

# Generate the contour plot 
contour(beta_seq, nu_seq, z, 
        levels = seq(0.01, 0.99, length.out = 25),
        drawlabels = FALSE,
        xlab = expression(beta), 
        ylab = expression(ln(sigma^2)),
        main = "Contour Plot of Joint Posterior Density")

##########################################
