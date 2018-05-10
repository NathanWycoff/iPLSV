#!/usr/bin/Rscript
#  R/checks/num_max_grads.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.09.2018

## Test analytic gradients of the incomplete posterior.

# Calculate using finite differences.

## Check that fixing values doe what we want.
#source('R/num_max.R')
#source('R/generate.R')
#source('R/lib.R')
require(iplsv)
source('R/num_max.R')

set.seed(1234)

K <- 3
V <- 4
M <- 2
N.mu <- 2
P <- 2
eta <- 1
gamma <- 0.1 * K
beta <- 0.1 * M

ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

# Convert docs from TF form to a list of vecs
docs_list <- lapply(1:M, function(i) unlist(sapply(1:V, function(j) rep(j, ret$docs[i,j]))))
docs <- docs_list

h <- 1e-6
nP <- matrix(NA, nrow = K, ncol = P)
for (k in 1:K) {
    for (p in 1:P) {
        f <- nclp(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        ret$PSI[k,p] <- ret$PSI[k,p] + h
        fp <- nclp(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        nP[k,p] <- (fp - f) / h
    }
}

nT <- matrix(NA, nrow = M, ncol = P)
for (m in 1:M) {
    for (p in 1:P) {
        f <- nclp(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        ret$THETA[m,p] <- ret$THETA[m,p] + h
        fp <- nclp(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        nT[m,p] <- (fp - f) / h
    }
}

nH <- matrix(NA, nrow = K, ncol = V)
for (k in 1:K) {
    for (v in 1:V) {
        f <- nclp(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        ret$PHI[k,v] <- ret$PHI[k,v] + h
        fp <- nclp(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        nH[k,v] <- (fp - f) / h
    }
}

g_nilp(ret$PHI, ret$THETA, ret$PSI, docs, eta, gamma, beta, soft_PHI = FALSE)
print(list(nP, nT, nH))

## PHi grad one term at a time
phi_w <- c(0.7, 0.3)
rho_i <- c(0.2, 0.8)

ll <- function(phi_w) {
    PI <- sum(phi_w * rho_i)
    log(PI)
}

dll <- function(phi_w) {
    PI <- sum(phi_w * rho_i)
    rho_i[1] / PI[1]
}

phi_w_m <- phi_w
phi_w_m[1] <- phi_w_m[1] + h
(ll(phi_w_m) - ll(phi_w)) / h
dll(phi_w)
