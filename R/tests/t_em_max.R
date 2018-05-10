#!/usr/bin/Rscript
#  R/tests/t_em_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

require(iplsv)
require(LdaFuncs)
require(ggplot2)
source("R/lib.R")

set.seed(123)

K <- 3
V <- 4000
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

#Test using smart init
set.seed(123)
system.time(fit <- em_plsv(docs, K, V, P, eta, gamma, beta, 
                          make_plot = FALSE, THETA_init = 'smart', 
                          PSI_init = 'smart', PHI_init = 'smart', 
                          THETA_fix = list(), PSI_fix = list(),
                          verbose = TRUE, thresh = 5e-2,
                          max_iters = 1e3, lik_grad = 'Cpp'))

comp_scat2d(ret$THETA, ret$PSI)
quartz()
comp_scat2d(fit$par$THETA, fit$par$PSI)

# Test that we don't move away from the truth
set.seed(123)
system.time(fit <- em_plsv(docs, K, V, P, eta, gamma, beta, 
                          make_plot = FALSE, THETA_init = ret$THETA, 
                          PSI_init = ret$PSI, PHI_init = ret$PHI, 
                          THETA_fix = list(), PSI_fix = list(),
                          verbose = TRUE, thresh = 5e-2,
                          max_iters = 1e3, lik_grad = 'Cpp'))

comp_scat2d(ret$THETA, ret$PSI)
quartz()
comp_scat2d(fit$par$THETA, fit$par$PSI)
