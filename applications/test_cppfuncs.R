#!/usr/bin/Rscript
#  applications/test_cppfuncs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.08.2018

## Test cpp implementations for speed and accuracy.
require(Rcpp)
require(iplsv)
require(microbenchmark)
source('R/lib.R')
source('R/em_max.R')

sourceCpp("src/exp_nlpost.cpp")

K <- 3
V <- 4
M <- 20
N.mu <- 30
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M

ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

# Convert docs from TF form to a list of vecs
docs_list <- lapply(1:M, function(i) unlist(sapply(1:V, function(j) rep(j, ret$docs[i,j]))))
docs <- docs_list

est <- ret

##### E STEP
Ns <- sapply(docs, length)
Z_exp <- lapply(1:M, function(i) matrix(NA, nrow = Ns[i], ncol = K))

# Calculate topic membership for each doc
RHO <- matrix(NA, nrow = M, ncol = K)
for (i in 1:M) {
    dists <- sapply(1:K, function(k) sum((est$THETA[i,] - est$PSI[k,])^2))
    RHO[i,] <- softmax(-0.5 * dists)
}

# Propogate to E step
for (i in 1:M) {
    for (j in 1:Ns[i]) {
        w <- docs[[i]][j]
        for (k in 1:K) {
            Z_exp[[i]][j,k] <- RHO[i,k] * est$PHI[k,w]
        }
    }
    Z_exp[[i]] <- Z_exp[[i]] / rowSums(Z_exp[[i]])
}


microbenchmark(exp_nlpostC(Z_exp, ret$PHI, ret$THETA, ret$PSI, docs, Ns, eta, gamma, beta))
microbenchmark(exp_nlpost(Z_exp, ret$PHI, ret$THETA, ret$PSI, docs, eta, gamma, beta))
