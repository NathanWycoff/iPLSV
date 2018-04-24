#!/usr/bin/Rscript
#  R/checks/inits.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.24.2018

## Convince myself that inits work.

set.seed(1234)

# Generate some data
K <- 3
V <- 4
M <- 20
N.mu <- 3000
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M

ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)


## Do a first estimate with random inits
fit <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta)

## Do another, but init on the old one. Should be instantaneous optim and no change to params.
fit_i <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, THETA_init = fit$par$THETA, 
                       PSI_init = fit$par$PSI, PHI_init = fit$par$PHI)
