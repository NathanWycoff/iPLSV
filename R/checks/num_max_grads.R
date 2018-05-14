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
source('R/lib.R')

set.seed(1234)

K <- 3
V <- 4
M <- 20
N.mu <- 20
P <- 2
eta <- 1
gamma <- 0.1 * K
beta <- 0.1 * M

ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

h <- 1e-6
nP <- matrix(NA, nrow = K, ncol = P)
for (k in 1:K) {
    for (p in 1:P) {
        f <- nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        ret$PSI[k,p] <- ret$PSI[k,p] + h
        fp <- nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        nP[k,p] <- (fp - f) / h
    }
}

nT <- matrix(NA, nrow = M, ncol = P)
for (m in 1:M) {
    for (p in 1:P) {
        f <- nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        ret$THETA[m,p] <- ret$THETA[m,p] + h
        fp <- nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        nT[m,p] <- (fp - f) / h
    }
}

nH <- matrix(NA, nrow = K, ncol = V)
for (k in 1:K) {
    for (v in 1:V) {
        f <- nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        ret$PHI[k,v] <- ret$PHI[k,v] + h
        fp <- nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
        nH[k,v] <- (fp - f) / h
    }
}

# Numerical gradient between PHI and GAMMA, the latent space PHI.
GAMM <- inv_softmax(ret$PHI)
nPG <- matrix(NA, nrow = K, ncol = V-1)
for (k in 1:K) {
    for (v in 1:(V-1)) {
        phi <- t(apply(cbind(GAMM, 0), 1, softmax))[1,1]
        GAMM[k,v] <- GAMM[k,v] + h
        phip <- t(apply(cbind(GAMM, 0), 1, softmax))[1,1]
        nPG[k,v] <- (phip - phi) / h
    }
}

# try to predict it 

nR <- matrix(NA, nrow = K, ncol = V-1)
ret$PHI <- inv_softmax(ret$PHI)
for (k in 1:K) {
    for (v in 1:(V-1)) {
        f <- nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE)
        ret$PHI[k,v] <- ret$PHI[k,v] + h
        fp <- nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE)
        nR[k,v] <- (fp - f) / h
    }
}

g_nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = FALSE)
g_nlip(ret$PHI, ret$THETA, ret$PSI, ret$docs, eta, gamma, beta, soft_PHI = TRUE)
print(list(nP, nT, nH))
