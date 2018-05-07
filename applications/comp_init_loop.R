#!/usr/bin/Rscript
#  applications/comp_init_loop.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

## Compare between multiple initialization methods.

require(iplsv)
require(LdaFuncs)

## Simulation Params
data_sets <- 15#How many datasets to simulate?
rand_inits <- 1#For each dataset, how many random inits?
set.seed(123)

## Data params
K <- 3
V <- 4000
M <- 20
N.mu <- 300
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M


ldapca_times <- rep(NA, data_sets)
ldapca_costs <- rep(NA, data_sets)
rand_times <- rep(NA, data_sets)
rand_costs <- rep(NA, data_sets)
for (d in 1:data_sets) {
    print(d)

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
                              verbose = FALSE, thresh = 5e-2,
                              max_iters = 1e3))
    ldapca_times[d] <- tt[1]
    ldapca_costs[d] <- nclp(fit$par$PHI, fit$par$THETA, fit$par$PSI, ret$docs, eta, gamma, beta)


    # Do a bunch of random inits
    d_times <- rep(NA, rand_inits)
    d_costs <- rep(NA, rand_inits)
    for (i in 1:rand_inits) {
        tt <- system.time(fit <- em_plsvn(docs, K, V, P, eta, gamma, beta, 
                                  make_plot = FALSE, 
                                  THETA_fix = list(), PSI_fix = list(),
                                  verbose = FALSE, thresh = 5e-2,
                                  max_iters = 1e3))
        d_times[i] <- tt[1]
        d_costs[i] <- nclp(fit$par$PHI, fit$par$THETA, fit$par$PSI, ret$docs, eta, gamma, beta)
    }

    rand_times[d] <- mean(d_times)
    rand_costs[d] <- min(d_costs)
}

save(rand_times, rand_costs, ldapca_time, ldapca_cost, file = 'data/comp_init.RData')
