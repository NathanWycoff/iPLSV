#!/usr/bin/Rscript
#  R/tests/em_vs_nilp.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.09.2018

## Compare using the EM algo to numerical maximization of the negative incomplete log posterior.

require(iplsv)
require(LdaFuncs)
require(ggplot2)

set.seed(1234)

K <- 3
V <- 4
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

THETA_init <- matrix(rnorm(M*P), nrow = M, ncol = P)
PSI_init <- matrix(rnorm(K*P), nrow = K, ncol = P)
PHI_init <- matrix(rgamma(K*V, 1, 1), nrow = K, ncol = V)
PHI_init <- PHI_init / rowSums(PHI_init)

set.seed(123)
fit_em <- em_plsv(docs, K, V, P, eta, gamma, beta, lik_grad = 'Cpp', thresh = 1e-3,
                  THETA_init = THETA_init,
                  PSI_init = PSI_init,
                  PHI_init = PHI_init, verbose = TRUE)

set.seed(123)
fit_num <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, 
                  #THETA_init = THETA_init,
                  #PSI_init = PSI_init,
                  #PHI_init = PHI_init)
                  THETA_init = fit_em$par$THETA,
                  PSI_init = fit_em$par$PSI,
                  PHI_init = fit_em$par$PHI)

comp_scat2d(fit_em$par$THETA, fit_em$par$PSI)
quartz()
comp_scat2d(fit_num$par$THETA, fit_num$par$PSI)
comp_scat2d(THETA_init, PSI_init)
