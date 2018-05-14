#!/usr/bin/Rscript
#  R/checks/check_new_lik.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.13.2018

#Test out the num_max function.
set.seed(1234)
require(iplsv)
source("R/num_max.R")
source("R/generate.R")
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

## LIkelihood should differ only by a constant
o1 <- nlip(ret$PHI, ret$THETA, ret$PSI, docs, eta, gamma, beta, soft_PHI = FALSE, new = FALSE)
n1 <- nlip(ret$PHI, ret$THETA, ret$PSI, docs_list, eta, gamma, beta, soft_PHI = FALSE, new = TRUE)

new_PSI <- matrix(rnorm(K*P), ncol = P)
o2 <- nlip(ret$PHI, ret$THETA, new_PSI, docs, eta, gamma, beta, soft_PHI = FALSE, new = FALSE)
n2 <- nlip(ret$PHI, ret$THETA, new_PSI, docs_list, eta, gamma, beta, soft_PHI = FALSE, new = TRUE)


## We should get the same solution with similar inits
fit1 <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, PHI_init = ret$PHI, THETA_init = ret$THETA, PSI_init = ret$PSI)
