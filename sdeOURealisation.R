# create N SDEs with gamma = 2, mu = -1
# and use MLE to find estimates of gamma and mu

N <- 500
gamma <- 7
mu <- -1
t <- seq(from = 0, to = 10, length = 10000)
sdeList <- vector(mode = "list", length = N)
sdeList <- lapply(sdeList, function(x) {
                          ouProcessEMApproximation(gamma = gamma, mu = mu, tau = t)$y
        })

sdeMLE <- t(sapply(sdeList, function(y) {
                         mleTmp <- mleOU(X=y, t=t, gamma0 = 1, mu0 = 1, verbose = FALSE)
                         c(mleTmp$gamma, mleTmp$mu)
        }))


mleMean <- colSums(sdeMLE) / (N - 1)



mleOU <- function(X, t, gamma0, mu0, max.iter=1000, diff.tol = 1e-8) {

