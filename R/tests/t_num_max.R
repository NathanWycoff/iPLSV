#!/usr/bin/Rscript
#  R/tests/t_num_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.22.2018

#Test out the num_max function.
set.seed(1234)
require(iplsv)
source("R/num_max.R")
require(ggplot2)

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
fit2 <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, PHI_init = ret$PHI, THETA_init = ret$THETA, PSI_init = ret$PSI)

comp_scat2d(fit1$par$THETA, fit1$par$PSI)
quartz()
comp_scat2d(fit2$par$THETA, fit2$par$PSI)

comp_scat2d(ret$THETA, ret$PSI)
