euler_maruyama(f, g, x0, t0, t1, dt){
    # f: drift function, must take x and t as arguments
    # g: diffusion function, must take x and t as arguments
    # x0: initial condition, matrix of size (B, d)
    # t0: initial time, scalar
    # t1: final time, scalar
    # dt: time step
    # returns: list of times and positions (also a list, with component having siez (B, d))
    n <- round((t1 - t0) / dt)
    x <- numeric(n)
    t <- numeric(n)
    x <- list()
    x[[1]] <- x0
    t[1] <- t0
    for (i in 2:n){
        t[i] <- t[i - 1] + dt
        x[[i]] <- x[[i - 1]] + f(x[[i - 1]], t[i - 1]) * dt + g(x[[i - 1]], t[i - 1]) * rnorm(1) * sqrt(dt)
    }
    return(list(t = t, x = x))
}



