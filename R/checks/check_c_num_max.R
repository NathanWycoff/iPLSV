#!/usr/bin/Rscript
#  R/checks/check_c_num_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.14.2018

### Test that Rcpp's implementation of the num max funcs mathces the R implementation.
require(Rcpp)
require(RcppArmadillo)
sourceCpp('./src/nlip.cpp')
source("R/num_max.R")
require(iplsv)
library(microbenchmark)

set.seed(1234)

K <- 3
V <- 4
M <- 20
N.mu <- 20
P <- 2
eta <- 100
gamma <- 0.1 * K
beta <- 0.1 * M

# Test the likelihood
ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

PHI_n <- t(apply(ret$PHI, 1, inv_softmax))[,-V]
microbenchmark(nlip(PHI_n, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE),
nlipC(PHI_n, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta))

# Test the gradient
PHI_n <- t(apply(ret$PHI, 1, inv_softmax))[,-V]
Ns <- sapply(ret$docs, length)

microbenchmark(g_nlip(PHI_n, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE),
g_nlipC(PHI_n, ret$THETA, ret$PSI, Ns, ret$docs, eta, gamma, beta))
