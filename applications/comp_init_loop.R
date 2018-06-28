#!/usr/bin/Rscript
#  applications/comp_init_loop.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

## Compare between multiple initialization methods.

require(iplsv)
source('R/smart_init.R')
require(LdaFuncs)
require(parallel)

## Simulation Params
data_sets <- 30#How many datasets to simulate?
rand_inits <- 1#For each dataset, how many random inits?
set.seed(123)
seeds <- sample(1:1e6, data_sets)

## Data params
K <- 3
M <- 20
N.mu <- 300
V <- ceiling(44 * (M*N.mu)^0.5)# Heap's law
P <- 2
eta <- 2
gamma <- 0.1 * K
beta <- 0.1 * M

# Init parallelization 
no_cores <- detectCores() - 1
no_cores <- 4
cl <- makeCluster(no_cores, 'FORK')

ldapca_times <- rep(NA, data_sets)
ldapca_costs <- rep(NA, data_sets)
rand_times <- rep(NA, data_sets)
rand_costs <- rep(NA, data_sets)
ret <- parLapply(cl, 1:data_sets, function(d)  {
                     set.seed(seeds[d])
                     ret <- gen_plsv(K, V, M, N.mu, P, eta, gamma, beta)

                     # Convert docs from TF form to a list of vecs
                     docs <- ret$docs

                     # Do a smart initialization
                     init <- smart_init(docs, alpha = 1, eta, K, V)
                     tt <- system.time(fit <- num_post_plsv(docs, K, V, P, eta, gamma, beta, 
                                     make_plot = FALSE, 
                                     PHI_init = ret$PHI, THETA_init = ret$THETA, PSI_init = ret$PSI))

                     ldatime <- tt[1]
                     ldacost <- -fit$logpost


                     # Do a bunch of random inits
                     d_times <- rep(NA, rand_inits)
                     d_costs <- rep(NA, rand_inits)
                     for (i in 1:rand_inits) {
                         tt <- system.time(fit <- num_post_plsv(docs, K, V, P, eta, gamma, beta, 
                                         make_plot = FALSE))

                         d_times[i] <- tt[1]
                         d_costs[i] <- -fit$logpost
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

ggdf <- as.data.frame(cbind(c(rep('Smart', data_sets), rep('Random', data_sets)), c(ldapca_times, rand_times), c(ldapca_costs, rand_costs)))
colnames(ggdf) <- c('Type', 'Time', 'Cost')

save(rand_times, rand_costs, ldapca_times, ldapca_costs, ggdf, file = 'data/comp_init.RData')
