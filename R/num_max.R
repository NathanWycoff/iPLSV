#!/usr/bin/Rscript
#  R/num_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.20.2018

## Numerically maximize the posterior with latent vars integrated out.

# Negative Complete Log Posterior
nclp <- function(PHI_n, THETA, PSI, docs, eta, gamma, beta) {
    #Transform back to the simplex.
    PHI_n <- cbind(PHI_n, 0)
    PHI <- t(apply(PHI_n, 1, softmax))

    ## Generate parameters
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

    ## Transform Parameters
    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    # Prob of each word in each doc
    PI <- RHO %*% PHI

    ## Generate Data
    for (i in 1:M) {
        ll <- ll + dmultinom(docs[i,], prob = PI[i,], log = TRUE)
    }

    return(-ll)
}

#' Numerically Maximize complete Posterior.
#'
#' Numerically maximize the complete posterior of the PLSV model.
#'
#' @param docs A term frequency matrix, that is, one with a row for each document, a column for each vocab word, and integer entries indicating the occurence of a word in a doc.
#' @param K The number of topics, an integer scalar.
#' @param V The number of unique words, an integer scalar.
#' @param P The dimensionality of the embedding space, an integer, usually 2.
#' @param eta The exchangible dirichlet prior on words in a topic.
#' @param beta The precision for topic locations, a positive scalar.
#' @param gama The precision for document locations, a positive scalar.
#' @param make_plot A boolean, if TRUE, will make a ggplot visualization of the topics and documents, with topics in red.
#' @param THETA_fix A list of lists, used to fix rows of THETA to a given value. Each sublist has two elements: 'ind' and 'val'. 'ind' Indicates the row, 1-index, of THETA to fix, and 'val', a real valued P-vector, indicates the value to fix it to.
#' @return A list containing ests, a list with PHI, the topic by document matrix, THETA, the document locations in P-D space, and PSI, the topic locations in P-D space.
#' @export
num_post_plsv <- function(docs, K, V, P, eta, gamma, beta, 
                          make_plot = FALSE, THETA_fix = list(), PSI_fix = list()) {

    M <- nrow(docs)

    THETA_inds <- sapply(THETA_fix, function(i) i$ind)
    ## Maximize the log posterior with respect to the doc and topic locations.

    # A vector representation of our three matrices, for use in numerical optimization
    # We drop the rows which are fixed here
    mat3par <- function(PHI_n, PSI, THETA) {
        #Check if there's any fixed vals.
        if (length(THETA_inds) > 0) {
            THETA_pass <- as.numeric(THETA[-THETA_inds,])
        } else {
            THETA_pass <- as.numeric(THETA)
        }
        par <- c(as.numeric(PHI_n), 
                 #as.numeric(PSI[-PSI_inds,]), 
                 as.numeric(PSI), 
                 THETA_pass)
        return(par)
    }

    # Return back to the matrix repr, for use in numerical optimization
    par3mat <- function(par, K, V, P) {
        #Convert from vectors to matrices
        to1 <- K*(V-1)
        to2 <- K*P
        PHI_n <- par[1:to1]
        PSI <- par[(to1+1):(to1+to2)]
        THETA <- par[(to1+to2+1):(length(par))]

        #Recover the matrix form
        PHI_n <- matrix(PHI_n, ncol = V-1)
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

        return(list(PHI_n = PHI_n, THETA = THETA, PSI = PSI))
    }

    joint_nclp_wrap <- function(par) {
        ret <- par3mat(par, K, V, P)
        nclp(ret$PHI_n, ret$THETA, ret$PSI, docs, eta, gamma, beta)
    }


    #Initialize at the truth
    #PSI_start <- PSI
    #THETA_start <- THETA
    #PHI_n_start <- t(apply(PHI, 1, inv_softmax))[,-V]

    #Initialize randomly
    PSI_start <- matrix(rnorm(K*P), ncol = P)
    THETA_start <- matrix(rnorm(M*P), ncol = P)
    PHI_n_start <- matrix(rnorm(K*(V-1)), ncol = V-1)

    #TODO: Init intelligently.

    # Do the actual optimization
    fit <- optim(mat3par(PHI_n_start, PSI_start, THETA_start), joint_nclp_wrap, 
                 method = 'BFGS')
    ests <- par3mat(fit$par, K, V, P)

    #Get PHI from PHI_n, with a 0 fixed as the last word for identifiability.
    ests$PHI <- t(apply(cbind(ests$PHI_n, 0), 1, softmax))
    ests$PHI_n <- NULL

    if (make_plot == TRUE) {
        comp_scat2d(ests$THETA, ests$PSI)
    }

    return(list(par = ests, logpost = -fit$value))
}
