# create Euler-Maruyama approximation for Ornstein-Uhlenbeck

# first, consider OU process with
# dX = \alpha X dt + dW

ouProcessEMApproximation <- function(x0=0, gamma, mu=0, sigma=1, tau) {
        # approximate Orstein-Uhlenbeck process
        # dX = \gamma (\mu - X) dt + \sigma dW
        # by Euler-Maruyama approximation
        # tau is time discretisation in [t0, T]
        # t0 = \tau_0 < \tau_1 < ... < \tau_N = T
        N <- length(tau) - 1
        y0 <- x0 
        y <- numeric(length = (N + 1))
        deltaTau <- tau[2:(N+1)] - tau[1:N]
        deltaW <- rnorm(N, mean = 0, sd = sqrt(deltaTau))
        for(n in 1:N) {
                termOne <- gamma * (mu - y[n]) * deltaTau[n]
                termTwo <- sigma * deltaW[n]
                y[n+1] <- y[n] + termOne + termTwo
        }
        list(t = tau, y = y)
}

