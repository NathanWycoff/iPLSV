#!/usr/bin/Rscript
#  R/tests/t_em_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

require(Rcpp)
source('R/generate.R')
source('R/lib.R')
source('R/num_max.R')
source('R/em_max.R')

set.seed(123)

K <- 5
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


set.seed(123)
system.time(fit <- em_plsvn(docs, K, V, P, eta, gamma, beta, 
                          make_plot = FALSE, THETA_init = 'smart', 
                          PSI_init = 'smart', PHI_init = 'smart', 
                          THETA_fix = list(), PSI_fix = list(),
                          verbose = TRUE, thresh = 5e-2,
                          max_iters = 1e3))

sourceCpp("src/exp_nlpost.cpp")
set.seed(123)
system.time(fit <- em_plsvn(docs, K, V, P, eta, gamma, beta, 
                          make_plot = FALSE, THETA_init = 'smart', 
                          PSI_init = 'smart', PHI_init = 'smart', 
                          THETA_fix = list(), PSI_fix = list(),
                          verbose = TRUE, thresh = 5e-2,
                          max_iters = 1e3, lik_grad = 'Cpp'))

comp_scat2d(ret$THETA, ret$PSI)
quartz()
comp_scat2d(fit$par$THETA, fit$par$PSI)
