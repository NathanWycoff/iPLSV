#!/usr/bin/Rscript
#  r/checks/c_num_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.20.2018

## Some checks we did when writing num_max.R

# Negative Complete Log Posterior
nclp <- function(PHI_n, THETA, PSI, docs, eta, gamma, beta) {
    #Transform back to the simplex.
    PHI_n <- cbind(PHI_n, 0)
    PHI <- t(apply(PHI_n, 1, softmax))

    ## Generate parameters
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

    ## Transform Parameters
    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    PI <- RHO %*% PHI

    ## Generate Data
    for (i in 1:M) {
        ll <- ll + dmultinom(docs[i,], prob = PI[i,], log = TRUE)
    }

    return(-ll)
}

#Optimize WRT PSI only
PSI2par <- function(PSI) {
    as.numeric(PSI)
}

par2PSI <- function(par, P) {
    matrix(par, ncol = P)
}

PSI_nclp_wrap <- function(par) {
    nclp(PHI, THETA, par2PSI(par, P), docs, eta, gamma, beta)
}

ran_init <- matrix(rnorm(K*P), ncol = P)
fit <- optim(PSI2par(PSI), PSI_nclp_wrap, method = 'BFGS')
comp_scat2d(par2PSI(fit$par, 2), PSI)

#Optimize WRT THETA only
THETA2par <- function(THETA) {
    as.numeric(THETA)
}

par2THETA <- function(par, P) {
    matrix(par, ncol = P)
}

THETA_nclp_wrap <- function(par) {
    nclp(PHI, par2THETA(par, P), PSI, docs, eta, gamma, beta)
}

ran_init <- matrix(rnorm(M*P), ncol = P)
fit <- optim(THETA2par(THETA), THETA_nclp_wrap, method = 'BFGS')
comp_scat2d(par2THETA(fit$par, 2), THETA)

## Join Optimization
## Maximize the log posterior with respect to the doc and topic locations.
# A vector representation of our two matrices, for use in numerical optimization
mat3par <- function(PHI_n, PSI, THETA) {
    par <- c(as.numeric(PHI_n), as.numeric(PSI), as.numeric(THETA))
    return(par)
}

# Return back to the matrix repr, for use in numerical optimization
par3mat <- function(par, K, V, P) {
    #Convert from vectors to matrices
    to1 <- K*(V-1)
    to2 <- K*P
    PHI_n <- par[1:to1]
    PSI <- par[(to1+1):(to1+to2)]
    THETA <- par[(to1+to2+1):(length(par))]

    PHI_n <- matrix(PHI_n, ncol = V-1)
    PSI <- matrix(PSI, ncol = P)
    THETA <- matrix(THETA, ncol = P)

    return(list(PHI_n = PHI_n, THETA = THETA, PSI = PSI))
}

joint_nclp_wrap <- function(par) {
    ret <- par3mat(par, K, V, P)
    nclp(ret$PHI_n, ret$THETA, ret$PSI, docs, eta, gamma, beta)
}


#ran_init <- matrix(rnorm(M*P), ncol = P)
PSI_start <- PSI
THETA_start <- THETA
PHI_n_start <- t(apply(PHI, 1, inv_softmax))[,-V]
fit <- optim(mat3par(PHI_n_start, PSI_start, THETA_start), joint_nclp_wrap, method = 'BFGS')
ests <- par3mat(fit$par, K, V, P)

#Compare the two scatterplots
quartz()
comp_scat2d(ests$THETA, ests$PSI)
quartz()
comp_scat2d(THETA, PSI)

# Compare the two topic-word tables.
t(apply(cbind(ests$PHI_n, 0), 1, softmax))
PHI
