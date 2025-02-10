# Fractional Brownian motion----
 
## Code adapted from: https://github.com/732jhy/fractionalBM/blob/master/R/fbm_sim.R
## The original code is for 1D fBM, but we can extend it to 2D, and self-define initial location

## The following code is for 2D fBM with Hurst parameter H, and initial location (x0, y0)

## t0 Starting time (default = 0)
## t1 Ending time (default = 1)
## dt Time step size (default = 0.01)
## H Hurst parameter (0 < H < 1)
## x0 Starting x-coordinate (default = 0)
## y0 Starting y-coordinate (default = 0)
# rho Cross-correlation between x and y components (-1 <= rho <= 1)
## method Simulation method (currently only 'cholesky' supported)

## return a list containing:
## - times: Sequence of time points
## - x: X-coordinates of FBM path
## - y: Y-coordinates of FBM path


FBM2D <- function(t0=0, t1=1, dt=0.01, H=0.5, x0=0, y0=0, rho=0, method='cholesky') {
  
  # H validation
  if (H <= 0 || H >= 1) stop("H must be in (0,1)")
  if (abs(rho) > 1) stop("rho must be between [-1,1]")
  
  # Update time
  n <- round((t1 - t0)/dt)
  actual_t1 <- t0 + n*dt
  T <- actual_t1 - t0
  times <- seq(t0, actual_t1, length.out = n+1)
  
  # Autocovariance function for FGN
  gamma <- function(k, H) 0.5*(abs(k-1)^(2*H) - 2*abs(k)^(2*H) + abs(k+1)^(2*H))
  
  # covariance matrix
  construct_cov_matrix <- function(n, H, rho) {
    Gamma <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        Gamma[i,j] <- gamma(abs(i-j), H)
      }
    }
    
    # Create block matrix
    block_top <- cbind(Gamma, rho*Gamma)
    block_bottom <- cbind(rho*Gamma, Gamma)
    rbind(block_top, block_bottom)
  }
  
  # Generate correlated FGN using Cholesky decomposition
  generate_correlated_fgn <- function(n, H, rho) {
    # Create covariance matrix
    Sigma <- construct_cov_matrix(n, H, rho)
    
    # Cholesky decomposition
    L <- chol(Sigma)
    
    # Generate random normals
    Z <- matrix(rnorm(2*n), ncol=2*n)
    
    # Correlated samples
    X <- Z %*% L
    list(x = X[1,1:n], y = X[1,(n+1):(2*n)])
  }
  
  # Generate fractional Gaussian noise
  fgn <- generate_correlated_fgn(n, H, rho)
  
  # Scale FGN and create FBM
  scale_factor <- T^H / n^H
  fbm_x <- cumsum(c(x0, fgn$x * scale_factor))
  fbm_y <- cumsum(c(y0, fgn$y * scale_factor))
  
  # Trim to match time points
  list(
    times = times,
    x = fbm_x[1:(n+1)],
    y = fbm_y[1:(n+1)]
  )
}

## Example usage
# fbm <- FBM2D(t0=0, t1=1, dt=0.001, H=0.7, x0=0, y0=0, rho=0.5)
# plot(fbm$x, fbm$y, type='l', main='2D Fractional Brownian Motion')