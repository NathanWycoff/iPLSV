#!/usr/bin/Rscript
#  R/em_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.07.2018

source('R/generate.R')
source('R/lib.R')
source('R/num_max.R')

# The expected negative log posterior
# Z_exp should be a list of matrices
exp_nlpost <- function(Z_exp, PHI, THETA, PSI, docs, eta, gamma, beta) {
    Ns <- sapply(docs, length)

    ## Prior distribution, does not depend on latent params.
    ll <- 0
    # Topic by Word Matrix
    for (k in 1:K){
        ll <- ll + (eta - 1) * sum(log(PHI[k,]))
    }

    # Topic Locations
    for (k in 1:K) {
        for (p in 1:P) {
            ll <- ll + dnorm(PSI[k,p],0,1/sqrt(beta), log = TRUE)
        }
    }

    # Doc Locations
    for (i in 1:M) {
        for (p in 1:P) {
            ll <- ll + dnorm(THETA[i,p],0,1/sqrt(gamma), log = TRUE)
        }
    }

    # Expected Likelihood
    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    # Sum it up for each doc
    for (m in 1:M) {
        for (n in 1:Ns[m]) {
            for (k in 1:K) {
                w <- docs[[m]][n]
                ll <- ll + Z_exp[[m]][n,k] * log(RHO[m,k] * PHI[k, w])
            }
        }
    }

    return(-ll)
}

# Get gradients of the expected negative log posterior.
g_enlp <- function(Z_exp, PHI, THETA, PSI, docs, eta, gamma, beta) {
    Ns <- sapply(docs, length)
    # Prior Contribution
    grad_THETA <- gamma * THETA
    grad_PSI <- beta * PSI

    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    # Liklihood contribution
    for (i in 1:M) {
        for (j in 1:Ns[i]) {
            for (k in 1:K) {
                g <- (Z_exp[[i]][j,k] - RHO[i,k]) * (THETA[i,] - PSI[k,])
                grad_THETA[i,] <- grad_THETA[i,] + g
                grad_PSI[k,] <- grad_PSI[k,] - g
            }
        }
    }

    return(list(grad_THETA = grad_THETA, grad_PSI = grad_PSI))
}

##Initialize some things
##PSI_hat <- matrix(rnorm(K*P), ncol = P)
##THETA_hat <- matrix(rnorm(M*P), ncol = P)
#PSI_hat <- ret$PSI
#THETA_hat <- ret$THETA
#PHI_hat <- ret$PHI
#Ns <- sapply(docs, length)
#est <- list(PSI = PSI_hat, 
#             THETA = THETA_hat,
#             PHI = PHI_hat)
#
#diff <- Inf
#thresh <- 1e-2#An absolute thresh in log space implies a relative thresh on the posterior; 1e-2 gives about 1%.
#max_iters <- 1e3
#last_true <- Inf
#iter <- 1
#while (iter < max_iters && diff > thresh) {
#    ##### E STEP
#    Z_exp <- lapply(1:M, function(i) matrix(NA, nrow = Ns[i], ncol = K))
#
#    # Calculate topic membership for each doc
#    RHO <- matrix(NA, nrow = M, ncol = K)
#    for (i in 1:M) {
#        dists <- sapply(1:K, function(k) sum((est$THETA[i,] - est$PSI[k,])^2))
#        RHO[i,] <- softmax(-0.5 * dists)
#    }
#
#    # Propogate to E step
#    for (i in 1:M) {
#        for (j in 1:Ns[i]) {
#            w <- docs[[i]][j]
#            for (k in 1:K) {
#                Z_exp[[i]][j,k] <- RHO[i,k] * est$PHI[k,w]
#            }
#        }
#        Z_exp[[i]] <- Z_exp[[i]] / rowSums(Z_exp[[i]])
#    }
#
#    ##### M STEP
#    # Estimate PHI, closed form
#    est$PHI <- matrix(eta, nrow = K, ncol = V)
#    for (i in 1:M) {
#        for (j in 1:Ns[i]) {
#            w <- docs[[i]][j]
#            for (k in 1:K) {
#                est$PHI[k,w] <- est$PHI[k,w] + Z_exp[[i]][j,k]
#            }
#        }
#    }
#    est$PHI <- est$PHI / rowSums(est$PHI)
#
#    # Need a bunch of wrappers since optim works only with vectors.
#    mat2par <- function(PSI, THETA) {
#        par <- c(as.numeric(PSI), as.numeric(THETA))
#        return(par)
#    }
#
#    par2mat <- function(par, K, P) {
#        PSI <- matrix(par[1:(K*P)], ncol = P)
#        THETA <- matrix(par[(K*P+1):length(par)], ncol = P)
#
#        return(list(PSI = PSI, THETA = THETA))
#    }
#
#    costwrap <- function(par) {
#        pars <- par2mat(par, K, P)
#        exp_nlpost(Z_exp, est$PHI, pars$THETA, pars$PSI, docs, eta, gamma, beta)
#    }
#
#    gradwrap <- function(par) {
#        pars <- par2mat(par, K, P)
#        grads <- g_enlp(Z_exp, est$PHI, pars$THETA, pars$PSI, docs, eta, gamma, beta)
#        mat2par(grads$grad_PSI, grads$grad_THETA)
#    }
#
#    # Do the actual optimization!
#    fit <- optim(mat2par(est$PSI, est$THETA), costwrap, gradwrap, method = 'BFGS')
#    last_est <- est
#    n_est <- par2mat(fit$par, K, P)
#    n_est$PHI <- est$PHI
#    est <- n_est
#
#    #Check convergence based on how much THETA and PSI are changing.
#    diff <- max(abs(est$THETA - last_est$THETA))
#    cat("Norm of diff:") 
#    print(diff)
#    
#    # Check the true posterior value.
#    PHI_n<- t(apply(est$PHI, 1, inv_softmax))[,-V]
#    true_val <- nclp(PHI_n, est$THETA, est$PSI, ret$docs, eta, gamma, beta)
#    last_true <- true_val
#    cat("True nlog post:")
#    print(true_val)
#}

#### Need a bunch of wrappers since optim works only with vectors.
mat2par <- function(PSI, THETA, PSI_inds, THETA_inds) {
    #Check if there's any fixed vals.
    if (length(THETA_inds) > 0) {
        THETA_pass <- as.numeric(THETA[-THETA_inds,])
    } else {
        THETA_pass <- as.numeric(THETA)
    }
    if (length(PSI_inds) > 0) {
        PSI_pass <- as.numeric(PSI[-PSI_inds,])
    } else {
        PSI_pass <- as.numeric(PSI)
    }

    return(c(PSI_pass, THETA_pass))
}

par2mat <- function(par, K, P, PSI_inds, THETA_inds) {
    #Convert from vectors to matrices
    to1 <- (K - length(PSI_inds))*P

    PSI <- par[1:to1]
    THETA <- par[(to1+1):(length(par))]

    #Recover the matrix form
    PSI <- matrix(PSI, ncol = P)
    THETA <- matrix(THETA, ncol = P)

    #Recover the fixed rows.
    THETA <- t(sapply(1:M, function(i) {
                          if (i %in% THETA_inds) {
                              return(THETA_fix[[which(THETA_inds==i)]]$val)
                          } else {
                              ind <- i - sum(THETA_inds < i)
                              return(THETA[ind,])
                          }
}))
    PSI <- t(sapply(1:K, function(i) {
                        if (i %in% PSI_inds) {
                            return(PSI_fix[[which(PSI_inds==i)]]$val)
                        } else {
                            ind <- i - sum(PSI_inds < i)
                            return(PSI[ind,])
                        }
}))


    return(list(PSI = PSI, THETA = THETA))
}


#' EM posterior maximization for PLSV.
#'
#' Perform the EM algorithm described in Iwata et al to estimate a PSLV model.
#'
#' @param docs A term frequency matrix, that is, one with a row for each document, a column for each vocab word, and integer entries indicating the occurence of a word in a doc.
#' @param K The number of topics, an integer scalar.
#' @param V The number of unique words, an integer scalar.
#' @param P The dimensionality of the embedding space, an integer, usually 2.
#' @param eta The exchangible dirichlet prior on words in a topic.
#' @param beta The precision for topic locations, a positive scalar.
#' @param gama The precision for document locations, a positive scalar.
#' @param make_plot A boolean, if TRUE, will make a ggplot visualization of the topics and documents, with topics in red.
#' @param THETA_init A real matrix with as many rows as docs has and P many columns, giving an initial value for THETA.
#' @param PSI_init A real matrix with K many rows and P many columns, giving an initial value for PSI
#' @param PHI_init A matrix of K many V-1-simplex valued rows, giving the initial value for PHI
#' @param THETA_fix A list of lists, used to fix rows of THETA to a given value. Each sublist has two elements: 'ind' and 'val'. 'ind' Indicates the row, 1-index, of THETA to fix, and 'val', a real valued P-vector, indicates the value to fix it to.
#' @return A list containing ests, a list with PHI, the topic by document matrix, THETA, the document locations in P-D space, and PSI, the topic locations in P-D space.
#' @export
em_plsv <- function(docs, K, V, P, eta, gamma, beta, 
                          make_plot = FALSE, THETA_init = NULL, 
                          PSI_init = NULL, PHI_init = NULL, 
                          THETA_fix = list(), PSI_fix = list(),
                          verbose = FALSE, thresh = 1e-2,
                          max_iters = 1e3) {
    M <- length(docs)
    Ns <- sapply(docs, length)

    # Get the indices of fixed params.
    THETA_inds <- sapply(THETA_fix, function(i) i$ind)
    PSI_inds <- sapply(PSI_fix, function(i) i$ind)
    #TODO: This code doesn't work lmao
    #if (max(THETA_inds) > M || max(PSI_inds) > K) {
    #    stop("Check the indices of your fixed values: one of them is larger than the corresponding matrix")
    #}

    ## Maximize the log posterior with respect to the doc and topic locations.
    # Do random inits if none were provided
    if (is.null(PSI_init)) {
        PSI_init <- matrix(rnorm(K*P), ncol = P)
    }
    if (is.null(THETA_init)) {
        THETA_init <- matrix(rnorm(M*P), ncol = P)
    }
    if (is.null(PHI_init)) {
        PHI_init <- matrix(rgamma(K*V, 1, 1), ncol = V)
    } 

    #TODO: Default inits based on LDA

    est <- list(PSI = PSI_init, 
                 THETA = THETA_init,
                 PHI = PHI_init)

    diff <- Inf
    last_true <- Inf
    iter <- 0

    while (iter < max_iters && diff > thresh) {
        iter <- iter + 1
        ##### E STEP
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

        # Wrap cost and gradient to work with vectors.
        costwrap <- function(par) {
            pars <- par2mat(par, K, P, PSI_inds, THETA_inds)
            exp_nlpost(Z_exp, est$PHI, pars$THETA, pars$PSI, docs, eta, gamma, beta)
        }

        gradwrap <- function(par) {
            pars <- par2mat(par, K, P, PSI_inds, THETA_inds)
            grads <- g_enlp(Z_exp, est$PHI, pars$THETA, pars$PSI, docs, eta, gamma, beta)
            mat2par(grads$grad_PSI, grads$grad_THETA, PSI_inds, THETA_inds)
        }


        # Do the actual optimization!
        fit <- optim(mat2par(est$PSI, est$THETA, PSI_inds, THETA_inds), costwrap, gradwrap, method = 'BFGS')
        last_est <- est
        n_est <- par2mat(fit$par, K, P, PSI_inds, THETA_inds)
        n_est$PHI <- est$PHI
        est <- n_est

        if (fit$convergence) {
            warning("Convergence Failure in em_plsv M step")
        }

        #Check convergence based on how much THETA and PSI are changing.
        diff <- max(abs(est$THETA - last_est$THETA))
        if (verbose) {
            cat("Iter ")
            cat(iter)
            cat("; norm of diff:") 
            cat(diff)
            cat("\n")
        }
    }

    if (max_iters <= iter) {
        warning("max_iters reached in em_plsv")
    }

    if (make_plot == TRUE) {
        comp_scat2d(est$THETA, est$PSI)
    }

    return(list(par = est, elogpost = -fit$value))
}
