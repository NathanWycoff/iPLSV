#!/usr/bin/Rscript
#  R/generate.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.20.2018

#' Generate from the PLSV model
#'
#'
#' @param K The number of topics, an integer scalar.
#' @param V The number of unique words, an integer scalar.
#' @param M The number of documents, an integer scalar.
#' @param N.mu The average number of words per document.
#' @param P The dimensionality of the embedding space, an integer, usually 2.
#' @param eta The exchangible dirichlet prior on words in a topic.
#' @param beta The precision for topic locations, a positive scalar.
#' @param gama The precision for document locations, a positive scalar.
#' @return A list containing the generated data, in particular, PHI is the topic by word matrix (rows sum to 1), PSI is a topic location matrix, docs is a list of generated documents, each entry denotes index of vocab, THETA is the document locations, and Z is a list of matrices, with a 1 denoting which topic the word belonged to.
#' @export
gen_plsv <- function(K, V, M, N.mu, P, eta, gamma, beta) {
    if (min(eta) < 0 || gamma < 0 || beta < 0) {
        stop("Hyperparameters for this model must be positive")
    }
    ## Generate parameters
    # Topic by Word Matrix
    PHI <- matrix(NA, nrow = K, ncol = V)
    for (k in 1:K){
        PHI[k,] <- rgamma(V, 1, eta)
    }
    PHI <- PHI / rowSums(PHI)

    # Doc Locations
    THETA <- matrix(NA, nrow = M, ncol = P)
    for (i in 1:M ){
        for (p in 1:P) {
            THETA[i,p] <- rnorm(1,0,1/sqrt(gamma))
        }
    }

    # Topic Locations
    PSI <- matrix(NA, nrow = K, ncol = P)
    for (k in 1:K) {
        for (p in 1:P) {
            PSI[k,p] <- rnorm(1,0,1/sqrt(beta))
        }
    }

    Ns <- rpois(M, N.mu - 1) + 1

    ## Transform Parameters
    # Get probability of topic in each doc
    RHO <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:M) {
        dists <- sapply(1:K, function(k) sum((THETA[i,] - PSI[k,])^2))
        RHO[i,] <- softmax(-0.5 * dists)
    }

    PI <- RHO %*% PHI

    ## Generate Data
    docs <- matrix(NA, nrow = M, ncol = V)
    for (i in 1:M) {
        docs[i,] <- rmultinom(1, Ns[i], PI[i,])
    }

    # Convert docs from TF form to a list of vecs
    docs_list <- lapply(1:M, function(i) unlist(sapply(1:V, function(j) rep(j, ret$docs[i,j]))))

    return(list(PHI = PHI, THETA = THETA, PSI = PSI, docs = docs_list))
}
