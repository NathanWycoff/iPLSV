#!/usr/bin/Rscript
#  R/comp_init.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

## Compare between multiple initialization methods.
source("R/em_max.R")
source('R/generate.R')
source('R/lib.R')
source('R/num_max.R')

## Simulation Params
rand_inits <- 30

set.seed(123)

# Data params
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

# Do a smart initialization
set.seed(123)
tt <- system.time(fit <- em_plsvn(docs, K, V, P, eta, gamma, beta, 
                          make_plot = FALSE, THETA_init = 'smart', 
                          PSI_init = 'smart', PHI_init = 'smart', 
                          THETA_fix = list(), PSI_fix = list(),
                          verbose = TRUE, thresh = 5e-2,
                          max_iters = 1e3))
ldapca_time <- tt[1]
ldapca_cost <- nclp(t(apply(fit$par$PHI, 1, inv_softmax))[,-V], fit$par$THETA, fit$par$PSI, ret$docs, eta, gamma, beta)


# Do a bunch of random inits
rand_times <- rep(NA, rand_inits)
rand_costs <- rep(NA, rand_inits)
for (i in 1:rand_inits) {
    tt <- system.time(fit <- em_plsvn(docs, K, V, P, eta, gamma, beta, 
                              make_plot = FALSE, 
                              THETA_fix = list(), PSI_fix = list(),
                              verbose = FALSE, thresh = 5e-2,
                              max_iters = 1e3))
    rand_times[i] <- tt[1]
    rand_costs[i] <- nclp(t(apply(fit$par$PHI, 1, inv_softmax))[,-V], fit$par$THETA, fit$par$PSI, ret$docs, eta, gamma, beta)

    print(i)
}

#Plot the two data.
ggdf <- as.data.frame(cbind(1, rand_times, rand_costs))
pt <- ggplot(ggdf, aes(x = V1, y= rand_times)) +
    geom_boxplot() + 
    geom_hline(yintercept = ldapca_time)
pc <- ggplot(ggdf, aes(x = V1, y= rand_costs)) +
    geom_boxplot() + 
    geom_hline(yintercept = ldapca_cost)
pc

grid.arrange(pt, pc, nrow = 1)
