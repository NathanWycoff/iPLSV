#!/usr/bin/Rscript
#  R/em_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

set.seed(123)
source('R/generate.R')
source('R/lib.R')
source('R/num_max.R')

K <- 3
V <- 40
M <- 20
N.mu <- 300
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
#PSI_hat <- matrix(rnorm(K*P), ncol = P)
#THETA_hat <- matrix(rnorm(M*P), ncol = P)
PSI_hat <- ret$PSI
THETA_hat <- ret$THETA
PHI_hat <- ret$PHI
Ns <- sapply(docs, length)
est <- list(PSI = PSI_hat, 
             THETA = THETA_hat,
             PHI = PHI_hat)

diff <- Inf
thresh <- 1e-2#An absolute thresh in log space implies a relative thresh on the posterior; 1e-2 gives about 1%.
max_iters <- 1e3
last_cost <- Inf
last_true <- Inf
iter <- 1
while (iter < max_iters && diff > thresh) {
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
                Z_exp[[i]][j,k] <- RHO[i,k] * est$PHI[k,w]
            }
        }
        Z_exp[[i]] <- Z_exp[[i]] / rowSums(Z_exp[[i]])
    }

    ##### M STEP
    # Estimate PHI, closed form
    est$PHI <- matrix(eta, nrow = K, ncol = V)
    for (i in 1:M) {
        for (j in 1:Ns[i]) {
            w <- docs[[i]][j]
            for (k in 1:K) {
                est$PHI[k,w] <- est$PHI[k,w] + Z_exp[[i]][j,k]
            }
        }
    }
    est$PHI <- est$PHI / rowSums(est$PHI)

    # Need a bunch of wrappers since optim works only with vectors.
    mat2par <- function(PSI, THETA) {
        par <- c(as.numeric(PSI), as.numeric(THETA))
        return(par)
    }

    par2mat <- function(par, K, P) {
        PSI <- matrix(par[1:(K*P)], ncol = P)
        THETA <- matrix(par[(K*P+1):length(par)], ncol = P)

        return(list(PSI = PSI, THETA = THETA))
    }

    costwrap <- function(par) {
        pars <- par2mat(par, K, P)
        exp_nlpost(Z_exp, est$PHI, pars$THETA, pars$PSI, docs, eta, gamma, beta)
    }

    gradwrap <- function(par) {
        pars <- par2mat(par, K, P)
        grads <- g_enlp(Z_exp, est$PHI, pars$THETA, pars$PSI, docs, eta, gamma, beta)
        mat2par(grads$grad_PSI, grads$grad_THETA)
    }

    # Do the actual optimization!
    fit <- optim(mat2par(PSI_hat, THETA_hat), costwrap, gradwrap, method = 'BFGS')
    last_est <- est
    n_est <- par2mat(fit$par, K, P)
    n_est$PHI <- est$PHI
    est <- n_est

    # Check the true posterior value.
    PHI_n<- t(apply(est$PHI, 1, inv_softmax))[,-V]
    true_val <- nclp(PHI_n, est$THETA, est$PSI, ret$docs, eta, gamma, beta)
    diff <- abs(last_true - true_val)
    last_true <- true_val
    cat("True nlog post:")
    print(true_val)
}

comp_scat2d(ret$PSI, ret$THETA)
quartz()
comp_scat2d(est$PSI, est$THETA)
