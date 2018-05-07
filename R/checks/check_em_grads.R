#!/usr/bin/Rscript
#  R/checks/check_em_grads.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

set.seed(1234)
source('R/generate.R')
source('R/lib.R')
source('R/num_max.R')

K <- 3
V <- 4
M <- 20
N.mu <- 30
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M

ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

# Convert docs from TF form to a list of vecs
docs_list <- lapply(1:M, function(i) unlist(sapply(1:V, function(j) rep(j, ret$docs[i,j]))))
docs <- docs_list

#' The expected negative log posterior
exp_nlpost <- function(Z_exp, PHI, THETA, PSI, docs, eta, gamma, beta) {
    Ns <- sapply(docs, length)

    ## Prior distribution, does not depend on latent params.
    ll <- 0
    # Topic by Word Matrix
    for (k in 1:K){
        ll <- ll + (eta - 1) * sum(log(PHI[k,]))
    }

    # Topic Locations
    for (k in 1:K) {
        for (p in 1:P) {
            ll <- ll + dnorm(PSI[k,p],0,1/sqrt(beta), log = TRUE)
        }
    }

    # Doc Locations
    for (i in 1:M) {
        for (p in 1:P) {
            ll <- ll + dnorm(THETA[i,p],0,1/sqrt(gamma), log = TRUE)
        }
    }

    # Expected Likelihood
    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    # Sum it up for each doc
    for (m in 1:M) {
        for (n in 1:Ns[m]) {
            for (k in 1:K) {
                w <- docs[[m]][n]
                ll <- ll + Z_exp[[m]][n,k] * log(RHO[m,k] * PHI[k, w])
            }
        }
    }

    return(-ll)
}

# Get gradients of the expected negative log posterior.
g_enlp <- function(Z_exp, PHI, THETA, PSI, docs, eta, gamma, beta) {
    Ns <- sapply(docs, length)
    # Prior Contribution
    grad_THETA <- gamma * THETA
    grad_PSI <- beta * PSI

    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    # Liklihood contribution
    for (i in 1:M) {
        for (j in 1:Ns[i]) {
            for (k in 1:K) {
                g <- (Z_exp[[i]][j,k] - RHO[i,k]) * (THETA[i,] - PSI[k,])
                grad_THETA[i,] <- grad_THETA[i,] + g
                grad_PSI[k,] <- grad_PSI[k,] - g
            }
        }
    }

    return(list(grad_THETA = grad_THETA, grad_PSI = grad_PSI))
}

#Initialize some things
PSI_hat <- matrix(rnorm(K*P), ncol = P)
THETA_hat <- matrix(rnorm(M*P), ncol = P)
Ns <- sapply(docs, length)
est <- list(PSI = PSI_hat, 
             THETA = THETA_hat)

##### E STEP
Z_exp <- lapply(1:M, function(i) matrix(NA, nrow = Ns[i], ncol = K))

# Calculate topic membership for each doc
RHO <- matrix(NA, nrow = M, ncol = K)
for (i in 1:M) {
    dists <- sapply(1:K, function(k) sum((est$THETA[i,] - est$PSI[k,])^2))
    RHO[i,] <- softmax(-0.5 * dists)
}

# Propogate to E step
for (i in 1:M) {
    for (j in 1:Ns[i]) {
        w <- docs[[i]][j]
        for (k in 1:K) {
            Z_exp[[i]][j,k] <- RHO[i,k] * ret$PHI[k,w]
        }
    }
    Z_exp[[i]] <- Z_exp[[i]] / rowSums(Z_exp[[i]])
}

h <- 1e-6
nP <- matrix(NA, nrow = K, ncol = P)
for (k in 1:K) {
    for (p in 1:P) {
        f <- exp_nlpost(Z_exp, ret$PHI, est$THETA, est$PSI, docs, eta, gamma, beta)
        est$PSI[k,p] <- est$PSI[k,p] + h
        fp <- exp_nlpost(Z_exp, ret$PHI, est$THETA, est$PSI, docs, eta, gamma, beta)
        nP[k,p] <- (fp - f) / h
    }
}
nT <- matrix(NA, nrow = M, ncol = P)
for (m in 1:M) {
    for (p in 1:P) {
        f <- exp_nlpost(Z_exp, ret$PHI, est$THETA, est$PSI, docs, eta, gamma, beta)
        est$THETA[m,p] <- est$THETA[m,p] + h
        fp <- exp_nlpost(Z_exp, ret$PHI, est$THETA, est$PSI, docs, eta, gamma, beta)
        nT[m,p] <- (fp - f) / h
    }
}
g_enlp(Z_exp, ret$PHI, est$THETA, est$PSI, docs, eta, gamma, beta)
