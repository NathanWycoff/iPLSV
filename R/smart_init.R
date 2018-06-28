#!/usr/bin/Rscript
#  R/smart_init.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 06.28.2018

#' LDA based Initailization for PLSV
#'
#' First a Latent Dirichlet Allocation model is fit on the data. THETA, the document locations, is initialized from a PCA of the document-topic matrix given by the LDA fit, PSI, the topic locations, is a weighted average of the locations of all documents with weights given by the topic prevalence in each document. The topic by vocab word matrix is the same in both LDA and PSLV, and is used directly as an initial value.
#'
#' @param docs A list of integer vectors, each vector is a doc, each entry is the index of that word in the vocabulary.
#' @param alpha A vector of nonnegative reals; the prior on each topic in a document
#' @param eta A vector of nonnegative reals; the prior on each word in a topic.
#' @param K Number of topics
#' @param V Size of vocab
#' @return A list with components THETA_init, PSI_init, and PHI_init, matrices of appropriate sizes.
smart_init <- function(docs, alpha, eta, K, V) {
    # Fit an LDA model on the data.
    lda_fit <- wLDA(docs, alpha = alpha, eta = eta, K, V, iters = 1e3, weights = rep(1, V))

    THETA_init <- prcomp(lda_fit$GAMMA)$x[,1:2]
    norm_GAMM <- lda_fit$GAMMA / colSums(lda_fit$GAMMA)
    PSI_init <- t(norm_GAMM) %*% THETA_init
    PHI_init <- lda_fit$BETA

    return(list(THETA= THETA_init, PSI = PSI_init, PHI = PHI_init))
}
