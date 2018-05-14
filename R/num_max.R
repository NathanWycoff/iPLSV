#!/usr/bin/Rscript
#  R/num_max.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.20.2018

## Numerically maximize the posterior with latent vars integrated out.

#' Negative Log Incomplete Posterior
#'
#' The log posterior with the latent class assignments integrated out, evaluated at PHI, THETA, PSI as specified, given docs.
#'
#' @param PHI A K by V matrix with simplex valued rows giving the the probabilty of words (cols) in topics (rows).
#' @param THETA M by P real valued matrix, giving P-d locations of documents.
#' @param PSI K by P real valued matrix, giving P-d locations of topics
#' @param docs A term frequency matrix, that is, one with a row for each document, a column for each vocab word, and integer entries indicating the occurence of a word in a doc.
#' @param eta The exchangible dirichlet prior on words in a topic.
#' @param beta The precision for topic locations, a positive scalar.
#' @param gama The precision for document locations, a positive scalar.
#' @param soft_PHI A boolean, if TRUE, PHI was provided on the logodds scales.
#' @return A scalar, giving the negative log posterior density at the point.
#' @export
nlip <- function(PHI, THETA, PSI, docs, eta, gamma, beta, soft_PHI) {
    if (soft_PHI) {
        PHI <- t(apply(cbind(PHI, 0), 1, softmax))
    }

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
        doc <- docs[[i]]
        for (w in doc) {
            ll <- ll + log(PI[i, w])
        }
    }

    return(-ll)
}

#' The Gradient of the Incomplete Log Posterior
#' 
#' this boii accepts a list of docs for right now.
g_nilp <- function(PHI, THETA, PSI, docs, eta, gamma, beta, soft_PHI) {
    Ns <- sapply(docs, length)
    if (soft_PHI) {
        PHI <- t(apply(cbind(PHI, 0), 1, softmax))
    }

    # Prior Contribution
    grad_THETA <- -gamma * THETA
    grad_PSI <- -beta * PSI

    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    # Prob of each word in each doc
    PI <- RHO %*% PHI

    # Gradients for THETA
    for (i in 1:M) {
        for (j in 1:Ns[i]) {
            w <- docs[[i]][j]
            for (k in 1:K) {
                for (l in 1:K) {
                    grad_THETA[i,] <- grad_THETA[i,] + 
                        PHI[k, w] / PI[i, w] * RHO[i,k] * (as.numeric(l==k) - RHO[i,l]) * (PSI[l, ] - THETA[i,])
                }
            }
        }
    }

    #Gradients for PSI
    for (i in 1:M) {
        for (j in 1:Ns[i]) {
            w <- docs[[i]][j]
            for (k in 1:K) {
                for (l in 1:K) {
                    grad_PSI[k,] <- grad_PSI[k,] + 
                        PHI[l, w] / PI[i, w] * RHO[i,l] * (as.numeric(l==k) - RHO[i,k]) * (THETA[i,] - PSI[k,])
                }
            }
        }
    }

    #Gradients for PHI
    grad_PHI <- matrix(0, nrow = K, ncol = V)
    for (k in 1:K) {
        for (i in 1:M) {
            for (j in 1:Ns[i]) {
                w <- docs[[i]][j]
                print(RHO[i,k] / PI[i, w])
                grad_PHI[k, w] <- grad_PHI[k, w] + RHO[i,k] / PI[i, w]
            }
        }
    }

    return(list(grad_THETA = -grad_THETA, grad_PSI = -grad_PSI, grad_PHI = grad_PHI))
}

#' Numerically Maximize complete Posterior.
#'
#' Numerically maximize the complete posterior of the PLSV model.
#'
#' @param docs A list of integer vectors, each subvector representing a document in the form of the indices of its words in the vocab.
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
num_post_plsv <- function(docs, K, V, P, eta, gamma, beta, 
                          make_plot = FALSE, THETA_init = NULL, 
                          PSI_init = NULL, PHI_init = NULL, 
                          THETA_fix = list(), PSI_fix = list()) {


    M <- length(docs)

    # Get the indices of fixed params.
    THETA_inds <- sapply(THETA_fix, function(i) i$ind)
    PSI_inds <- sapply(PSI_fix, function(i) i$ind)
    #TODO: This code doesn't work lmao
    #if (max(THETA_inds) > M || max(PSI_inds) > K) {
    #    stop("Check the indices of your fixed values: one of them is larger than the corresponding matrix")
    #}

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
        if (length(PSI_inds) > 0) {
            PSI_pass <- as.numeric(PSI[-PSI_inds,])
        } else {
            PSI_pass <- as.numeric(PSI)
        }

        #Only otimize the fixed vals.
        par <- c(as.numeric(PHI_n), 
                 #as.numeric(PSI[-PSI_inds,]), 
                 PSI_pass, 
                 THETA_pass)
        return(par)
    }

    # Return back to the matrix repr, for use in numerical optimization
    par3mat <- function(par, K, V, P) {
        #Convert from vectors to matrices
        to1 <- K*(V-1)
        to2 <- (K - length(PSI_inds))*P

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
        PSI <- t(sapply(1:K, function(i) {
                    if (i %in% PSI_inds) {
                        return(PSI_fix[[which(PSI_inds==i)]]$val)
                    } else {
                        ind <- i - sum(PSI_inds < i)
                        return(PSI[ind,])
                    }
                 }))

        return(list(PHI_n = PHI_n, THETA = THETA, PSI = PSI))
    }

    joint_nlip_wrap <- function(par) {
        ret <- par3mat(par, K, V, P)
        nlip(ret$PHI_n, ret$THETA, ret$PSI, docs, eta, gamma, beta, soft_PHI = TRUE)
    }

    # Do random inits if none were provided
    if (is.null(PSI_init)) {
        PSI_init <- matrix(rnorm(K*P), ncol = P)
    }
    if (is.null(THETA_init)) {
        THETA_init <- matrix(rnorm(M*P), ncol = P)
    }
    if (is.null(PHI_init)) {
        PHI_n_init <- matrix(rnorm(K*(V-1)), ncol = V-1)
    } else {
        #Convert users simplex to real valued matrix for optim
        PHI_n_init <- t(apply(PHI_init, 1, inv_softmax))[,-V]
    }

    #TODO: Default inits based on LDA

    # Do the actual optimization
    fit <- optim(mat3par(PHI_n_init, PSI_init, THETA_init), joint_nlip_wrap, 
                 method = 'BFGS', control = list('maxit' = 1e3))
    ests <- par3mat(fit$par, K, V, P)

    if (fit$convergence) {
        warning("Convergence Failure in num_post_plsv")
    }

    #Get PHI from PHI_n, with a 0 fixed as the last word for identifiability.
    ests$PHI <- t(apply(cbind(ests$PHI_n, 0), 1, softmax))
    ests$PHI_n <- NULL

    if (make_plot == TRUE) {
        comp_scat2d(ests$THETA, ests$PSI)
    }

    return(list(par = ests, logpost = -fit$value))
}
