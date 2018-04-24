#!/usr/bin/Rscript
#  R/checks/fixing.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.22.2018

## Check that fixing values doe what we want.

set.seed(1234)

K <- 3
V <- 4
M <- 20
N.mu <- 3000
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M

ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

THETA_fix <- list(list(ind = 2, val = c(-1, 1)), list(ind = 7, val = c(-0.5, 0.5)))

## If we fix THETA to be exactly the truth, we ought to get the exact truth back from the optim routine.
THETA_fix <- lapply(1:M, function(i) list(ind = i, val = ret$THETA[i,]))

fit <- num_post_plsv(ret$docs, K, V, P, eta, gamma, beta, THETA_fix = THETA_fix)

fit$par$THETA
ret$THETA
