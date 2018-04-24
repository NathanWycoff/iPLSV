#!/usr/bin/Rscript
#  R/generate.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.20.2018

## Generate PLSV variates.

set.seed(1234)

K <- 3
V <- 4
M <- 20
N.mu <- 3000
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M

gen_plsv <- function(K, V, M, N.mu, P, eta, gamma, beta) {
    if (min(eta) < 0 || gamma < 0 || beta < 0) {
        stop("Hyperparameters for this model must be positive")
    }
    ## Generate parameters
    # Topic by Word Matrix
    PHI <- matrix(NA, nrow = K, ncol = V)
    for (k in 1:K){
        PHI[k,] <- rgamma(V, 1, eta)
    }
    PHI <- PHI / rowSums(PHI)

    # Doc Locations
    THETA <- matrix(NA, nrow = M, ncol = P)
    for (i in 1:M ){
        for (p in 1:P) {
            THETA[i,p] <- rnorm(1,0,1/sqrt(gamma))
        }
    }

    # Topic Locations
    PSI <- matrix(NA, nrow = K, ncol = P)
    for (k in 1:K) {
        for (p in 1:P) {
            PSI[k,p] <- rnorm(1,0,1/sqrt(beta))
        }
    }

    Ns <- rpois(M, N.mu - 1) + 1

    ## Transform Parameters
    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    PI <- RHO %*% PHI

    ## Generate Data
    docs <- matrix(NA, nrow = M, ncol = V)
    for (i in 1:M) {
        docs[i,] <- rmultinom(1, Ns[i], PI[i,])
    }

    return(list(PHI = PHI, THETA = THETA, PSI = PSI, docs = docs))
}
