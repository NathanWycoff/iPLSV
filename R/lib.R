#!/usr/bin/Rscript
#  R/lib.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 04.20.2018

require(ggplot2)

#A numerically stable softmax function
softmax <- function(x) {
    exp(x - max(x) - log(sum(exp(x - max(x)))))
}

#An inverse softmax function which assumes that the last argument of the softmax was zero.
inv_softmax <- function(rho) {
    a <- log(rho)
    a <- a - log(rho[length(rho)])
    return(a)
}

#'Solve the orthogonal procrustes problem
#'
#' @param A A Matrix of the same dimension as B
#' @param B A Matrix of the same dimension of A of which a rotation will be computed which puts it closest to A.
#' @param comp_stress Should we compute and return the F norm of A - WB, where W is the rotation?
#' @return If comp_stress, a list with elements W, the rotation matrix, and stress, the scalar norm of the  difference, or just the matrix W if not.
orth_proc <- function(A, B, comp_stress = TRUE) {
    ret <- svd(A %*% t(B))
    W <- ret$u %*% t(ret$v)

    if (comp_stress) {
        stress <- norm(A - W %*% B)
        return(list(W = W, stress = stress))
    } else {
        return(W)
    }
}



#Compare two sets of 2D points using a scatterplot
#'
#' Compare two sets of 2D points by plotting a scatterplot, optionally rotating the data such that they are as similar as possible (useful when visualizing projections which are invariante to unitary transformations (rotations and reflections). A will be in black, B will be in red.
#' @param A A first set of points, should be of the same dimensions as B
#' @param B A first set of points, should be of the same dimensions as A 
#' @param rot If TRUE, solves the orthogonal procrustes problem to rotate B to match A.
#' @param scale_ax If TRUE, scales each axis so as to minimize error (performed after procrustes analysis if rot is also TRUE).
#' @return
comp_scat2d <- function(A, B, rot = FALSE, scale_ax = FALSE, obnoxious = FALSE) {
    # Do procrustes on the transpose because we are rotating each row.
    if (rot) {
        W <- t(orth_proc(t(A), t(B), comp_stress = FALSE))
        B <- B %*% W
        if (obnoxious) {
            cat('Unitary Transform:')
            print(W)
        }
    }
    if (scale_ax) {
        s1 <- A[,1] %*% B[,1] / (B[,1] %*% B[,1])
        s2 <- A[,2] %*% B[,2] / (B[,2] %*% B[,2])
        B <- B %*% diag(c(s1, s2))
        if (obnoxious) {
            cat('Scaling:')
            print(c(s1, s2))
        }
    }

    # ggplot really likes dataframes
    Adf <- as.data.frame(A)
    Adf$ID <- rownames(Adf)
    Bdf <- as.data.frame(B)
    Bdf$ID <- rownames(Bdf)

    #Plot the two data.
    p <- ggplot(Adf, aes(x= V1, y= V2)) +
        geom_text(aes(label=ID)) +
        geom_text(data = Bdf,
            mapping = aes(x = V1, y = V2, label = ID, color = 'red'))

    print(p)
}

