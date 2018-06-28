#!/usr/bin/Rscript
#  R/tests/t_num_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.22.2018

#Test out the num_max function.
set.seed(1234)
require(iplsv)
source("R/num_max.R")
require(ggplot2)

## Test recovery of true locations.
K <- 3
V <- 4
M <- 20
N.mu <- 300
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M

ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

fit1 <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta)
fit2 <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta,
                      PHI_init = ret$PHI, THETA_init = ret$THETA, PSI_init = ret$PSI)

comp_scat2d(fit1$par$THETA, fit1$par$PSI)
dev.new()
comp_scat2d(fit2$par$THETA, fit2$par$PSI)

# Test smart initialization.
init <- smart_init(docs, alpha, eta, K, V)

fit <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, 
                      PHI_init = init$PHI, THETA_init = init$THETA, PSI_init = init$PSI)

# Before
comp_scat2d(init$THETA, init$PSI)
dev.new()
# After
comp_scat2d(fit$par$THETA, fit$par$PSI)

# Look at a random init
init <- list()
init$PSI <- matrix(rnorm(K*P), ncol = P)
init$THETA <- matrix(rnorm(M*P), ncol = P)
init$PHI <- matrix(rgamma(K*V, 1, 1), ncol = V)
init$PHI <- init$PHI / rowSums(init$PHI)

fit <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, 
                      PHI_init = init$PHI, THETA_init = init$THETA, PSI_init = init$PSI)

# Before
comp_scat2d(init$THETA, init$PSI)
dev.new()
# After
comp_scat2d(fit$par$THETA, fit$par$PSI)
