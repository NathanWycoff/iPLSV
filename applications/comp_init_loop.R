#!/usr/bin/Rscript
#  applications/comp_init_loop.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

## Compare between multiple initialization methods.

require(iplsv)
require(LdaFuncs)
require(parallel)

## Simulation Params
data_sets <- 30#How many datasets to simulate?
rand_inits <- 1#For each dataset, how many random inits?
set.seed(123)
seeds <- sample(1:1e6, data_sets)

## Data params
K <- 3
V <- 4000
M <- 20
N.mu <- 300
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M

# Init parallelization 
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, 'FORK')

ldapca_times <- rep(NA, data_sets)
ldapca_costs <- rep(NA, data_sets)
rand_times <- rep(NA, data_sets)
rand_costs <- rep(NA, data_sets)
ret <- parLapply(cl, 1:data_sets, function(d)  {
                     set.seed(seeds[d])
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
                     ldatime <- tt[1]
                     ldacost <- nclp(fit$par$PHI, fit$par$THETA, fit$par$PSI, ret$docs, eta, gamma, beta)


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

                     randtime <- mean(d_times)
                     randcost <- min(d_costs)

                     return(c(ldatime, ldacost, randtime, randcost))
})

#Store results
ldapca_times <- sapply(ret, function(i) i[1])
ldapca_costs <- sapply(ret, function(i) i[2])
rand_times <- sapply(ret, function(i) i[3])
rand_costs <- sapply(ret, function(i) i[4])

stopCluster(cl)

save(rand_times, rand_costs, ldapca_time, ldapca_cost, file = 'data/comp_init.RData')
