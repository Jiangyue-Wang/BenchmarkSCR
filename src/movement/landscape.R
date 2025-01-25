library(MASS)
library(matrix)

rmvn <- function(n, mu = 0, V = matrix(1)) {
    p <- length(mu)
    if (any(is.na(match(dim(V), p)))) 
        stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

ar1land <- function(grids = expand.grid(1:100, 1:100), phi = 0.05) {
    # grids: a matrix of size (n, 2) containing the coordinates of the grid points
    # phi: inverse length scale the AR(1) 
    
    nn <- nrow(grids)
    distance <- as.matrix(dist(grids))
    X <- rmvn(1, rep(0, nn), exp(-phi * distance))
    return(X)
}


