# want to estimate the log-likelihood and MLE for OU parameters
# dX = \gamma (\mu - X) dt + dW

ornsteinUhlenbeckProcess <- function(X, gamma, mu) {
        f <- gamma * (mu - X)
        f
}

stratonovichIntegral <- function(X, fnX) {
        n <- length(X)
        midpointf <- 0.5 * (fnX[1:(n-1)] + fnX[2:n])
        differenceX <- X[2:n] - X[1:(n-1)]
        S1n <- sum(midpointf * differenceX)
        S1n
}

# approximate the log-likelihood and MLE \gamma and \mu
.jntFunction <- function(X, t, mu, gamma) {
        n <- length(X)
        if(missing(t)) t <- 0:(n - 1)
        stopifnot(length(X) == length(t))
        terms <- (mu - X[1:(n-1)])^2 * (t[2:n] - t[1:(n-1)])
        gamma^2 * sum(terms)
}

.s1nFunction <- function(X, mu, gamma) {
        n <- length(X)
        terms <- (2 * mu - X[1:(n-1)] - X[2:n]) * (X[2:n] - X[1:(n-1)])
        0.5 * gamma * sum(terms)
}

likelihoodOU <- function(X, t, mu, gamma) {
        stopifnot(length(t) == length(X))
        n <- length(X)
        tn <- t[n]
        S1n <- .s1nFunction(X, mu, gamma)
        Jnt <- .jntFunction(X, t, mu, gamma)
        Lnt1 <- gamma * (S1n - 0.5 * tn) - 0.5 * gamma^2 * Jnt
        Lnt1
}

mleGammeOU <- function(X, t, mu) {
        # estimate MLE of gamme in 
        # dX = \gamma (\mu - X)dt + dW
        # with discrete obs from t0 = 0 to tn = T
        stopifnot(length(t) == length(X), 
                    t[1] == 0)
        n <- length(X)
        tn <- t[n]
        termOne <- (2 * mu - X[1:(n-1)] - X[2:n]) * (X[2:n] - X[1:(n-1)])
        termTwo <- 2 * (mu - X[1:(n-1)])^2 * (t[2:n] - t[1:(n-1)])
        top <- sum(termOne) + tn
        bottom <- sum(termTwo)
        gammaHat <- top / bottom
        gammaHat
}

mleMuOU <- function(X, t, gamma) {
        # estimate MLE of gamme in 
        # dX = \gamma (\mu - X)dt + dW
        # with discrete obs from t0 = 0 to tn = T
        stopifnot(length(t) == length(X), 
                    t[1] == 0)
        n <- length(X)
        tn <- t[n]
        XT <- X[n]
        X0 <- X[1]
        top <- XT - X0 + gamma * sum(X[1:(n-1)] * (t[2:n] - t[1:(n-1)]))
        bottom <- gamma * tn
        muHat <- top / bottom
        muHat
}

mleOU <- function(X, t, gamma0, mu0, max.iter=1000, diff.tol = 1e-8, 
                  verbose = TRUE) {
        # estimate MLE of gamme in 
        # dX = \gamma (\mu - X)dt + dW
        # with discrete obs from t0 = 0 to tn = T
        # using iterative scheme to find gamma and mu simultaneously
        stopifnot(length(t) == length(X), 
                    t[1] == 0)
        n <- length(X)
        tn <- t[n]
        i <- 0
        gamma <- numeric(length = max.iter)
        gamma[1] <- gamma0
        mu <- numeric(length = max.iter)
        mu[0] <- mu0
        diff <- 1
        while(i < max.iter & diff > diff.tol) {
                i <- i + 1
                gammaHat <- mleGammeOU(X=X, t=t, mu=mu[i])
                muHat <- mleMuOU(X=X, t=t, gamma=gammaHat)
                gamma[i+1] <- gammaHat
                mu[i+1] <- muHat
                diffGamma <- (gamma[i+1] - gamma[i])^2
                diffMu <- (mu[i+1] - mu[i])^2
                diff <- diffGamma + diffMu
        }
        if (verbose == TRUE) {
                cat("Completed in ", i, "iterations, with diff = ", diff, "\n")
        }
        res <- list(gammaHat = gamma[i+1], muHat = mu[i+1])
        res
}


              
