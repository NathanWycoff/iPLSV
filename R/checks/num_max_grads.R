#!/usr/bin/Rscript
#  R/checks/num_max_grads.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.09.2018

## Test analytic gradients of the incomplete posterior.

# Calculate using finite differences.

## Check that fixing values doe what we want.
#source('R/num_max.R')
#source('R/generate.R')
#source('R/lib.R')
require(iplsv)
source('R/num_max.R')
source('R/lib.R')

set.seed(1234)

K <- 3
V <- 4
M <- 20
N.mu <- 20
P <- 2
eta <- 100
gamma <- 0.1 * K
beta <- 0.1 * M

# Test the grads
ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

PHI_prior_grad <- function(PHI, THETA, PSI, docs, eta, gamma, beta, soft_PHI) {
    if (soft_PHI) {
        PHI <- t(apply(cbind(PHI, 0), 1, softmax))
    }

    grad_PHI <- (eta - 1) / PHI

    if (soft_PHI) {
        grad_GAMM <- matrix(0, nrow = K, ncol = V-1)
        for (k in 1:K) {
            for (v in 1:(V-1)) {
                grad_GAMM[k,v] <- (eta - 1) * (1 - V * PHI[k,v])
            }
        }
        return(grad_GAMM)
    } else {
        return(grad_PHI)
    }
}

PHI_prior <- function(PHI, THETA, PSI, docs, eta, gamma, beta, soft_PHI) {
    if (soft_PHI) {
        PHI <- t(apply(cbind(PHI, 0), 1, softmax))
    }

    nP <- matrix(0, nrow = K, ncol = V) 
    pp <- 0

    for (k in 1:K){
        pp <- pp + (eta - 1) * sum(log(PHI[k,]))
    }

    return(pp)
}

PHI_prior_num <- function(PHI, THETA, PSI, docs, eta, gamma, beta, soft_PHI) {

    nH <- matrix(NA, nrow = K, ncol = V)
    h <- 1e-8
    for (k in 1:K) {
        for (v in 1:V) {
            f <- PHI_prior(PHI, THETA, PSI, docs, eta, gamma, beta, soft_PHI = soft_PHI)
            PHI[k,v] <- PHI[k,v] + h
            fp <- PHI_prior(PHI, THETA, PSI, docs, eta, gamma, beta, soft_PHI = soft_PHI)
            nH[k,v] <- (fp - f) / h
        }
    }

    return(nH)
}

# try to predict it 
PHI_prior_grad(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
PHI_prior_num(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)

PHI_prior_grad(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE)
PHI_prior_num(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE)

g_nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
num_nlip_grad(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)


PHI_n <- t(apply(ret$PHI, 1, inv_softmax))[,-V]
g_nlip(PHI_n, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE)
num_nlip_grad(PHI_n, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE)

# Test the whole optim 
#set.seed(1234)
#K <- 3
#V <- 4
#M <- 20
#N.mu <- 30
#P <- 2
#eta <- 2
#gamma <- 0.1 * K
#beta <- 0.1 * M

set.seed(1234)

K <- 3
V <- 20
M <- 20
N.mu <- 10
P <- 2
eta <- 10
gamma <- 0.1 * K
beta <- 0.1 * M

ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

# Compare auto grads to an grads
system.time(fit1 <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, PHI_init = ret$PHI, THETA_init = ret$THETA, PSI_init = ret$PSI))
system.time(fit2 <- num_post_plsv_new(ret$docs, K, V, P, eta, gamma, beta, PHI_init = ret$PHI, THETA_init = ret$THETA, PSI_init = ret$PSI))

comp_scat2d(fit1$par$THETA, fit1$par$PSI)
quartz()
comp_scat2d(fit2$par$THETA, fit2$par$PSI)


# Compare auto grads to num grads
system.time(fit1 <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, PHI_init = ret$PHI, THETA_init = ret$THETA, PSI_init = ret$PSI))
system.time(fit2 <- num_post_plsv_new2(ret$docs, K, V, P, eta, gamma, beta, PHI_init = ret$PHI, THETA_init = ret$THETA, PSI_init = ret$PSI))

comp_scat2d(fit1$par$THETA, fit1$par$PSI)
quartz()
comp_scat2d(fit2$par$THETA, fit2$par$PSI)
