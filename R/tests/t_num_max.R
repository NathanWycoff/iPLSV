#!/usr/bin/Rscript
#  R/tests/t_num_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.22.2018

#Test out the num_max function.
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

fit <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta)
comp_scat2d(fit$par$THETA, fit$par$PSI)
