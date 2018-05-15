/** src/nlip.cpp Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.14.2018  */

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Note: parameterized in terms of phi = precision = 1/variance = 1/sqrt(std)
double logdnorm(double x, double mu, double phi) {
    // TODO: Remove normalizing constant.
    return 0.5 * log(phi / (2.0 * M_PI)) - 0.5 * pow(x - mu, 2) * phi;
}

// Numerically stable softmax function.
rowvec softmaxC(rowvec x) {
    return exp(x - max(x) - log(sum(exp(x - max(x)))));

}

// Turn PHI, a log odds matrix, rowwise into a simplex valued matrix (will gain 1 column).
mat mat_softmax(mat PHI) {
    int K = PHI.n_rows;
    int V = PHI.n_cols + 1;

    // Put PHI on the simplex scale.
    mat PHI_s(K,V);
    mat PHI_0 = join_rows(PHI, zeros(K));
    for (int k = 0; k < K; k++) {
        PHI_s.row(k) = softmaxC(PHI_0.row(k));
    }

    return(PHI_s);
}

// Get probability of topic in each doc
mat make_RHO(mat THETA, mat PSI, mat PHI_s) {
    int K = PSI.n_rows;
    int M = THETA.n_rows;

    mat RHO(M ,K);
    rowvec dists(K);
    for (int m = 0; m < M; m++) {
        for (int k = 0; k < K; k++) {
            dists(k) = sum(pow(THETA.row(m) - PSI.row(k), 2));
        }
        RHO.row(m) = softmaxC(-0.5 * dists);
    }

    return(RHO);
}

// Note: This assumes PHI is on the logodds scale, and has V many columns, not V-1 many.
//
// [[Rcpp::export]]
double nlipC(mat PHI, mat THETA, mat PSI, List docs, double eta, double gamma, double beta) {
    int K = PHI.n_rows;
    int V = PHI.n_cols + 1;
    int M = THETA.n_rows;
    int P = THETA.n_cols;

    mat PHI_s = mat_softmax(PHI);
    //// Generate parameters
    double ll = 0;

    //// Topic by Word Matrix
    for (int k = 0; k < K; k++) {
        ll += (eta - 1) * sum(log(PHI_s.row(k)));
    }

    // Topic Locations
    for (int k = 0; k < K; k++) {
        for (int p = 0; p < P; p++) {
            ll += logdnorm(PSI(k,p), 0, beta);
        }
    }

    // Doc Locations
    for (int i = 0; i < M; i++) {
        for (int p = 0; p < P; p++) {
            ll += logdnorm(THETA(i,p), 0, gamma);
        }
    }

    //// Transform Parameters
    mat RHO = make_RHO(THETA, PSI, PHI_s);
    mat ETA = RHO * PHI_s;

    //// Generate Data
    uvec doc;
    for (int i = 0; i < M; i++) {
        doc = as<uvec>(docs[i]);
        for (int w : doc) {
            ll += log(ETA(i, w-1));
        }
    }

    return(-ll);
}

double scalar_invert(double x) { return(1/x); }

// [[Rcpp::export]]
List g_nlipC(mat PHI, mat THETA, mat PSI, uvec Ns, List docs, double eta, double gamma, 
        double beta) {

    int K = PHI.n_rows;
    int V = PHI.n_cols + 1;
    int M = THETA.n_rows;
    int P = THETA.n_cols;

    mat PHI_s = mat_softmax(PHI);

    // Prior Contribution
    mat grad_THETA(THETA);
    grad_THETA *= -gamma;

    mat grad_PSI(PSI);
    grad_PSI *= -beta;

    mat grad_PHI(K, V, fill::zeros);
    for (int k = 0; k < K; k++) {
        for (int v = 0; v < V; v++) {
            grad_PHI(k,v) = (eta - 1) / PHI_s(k,v);
        }
    }

    // Transform Parameters
    mat RHO = make_RHO(THETA, PSI, PHI_s);
    mat ETA = RHO * PHI_s;


    // Gradients for THETA
    int w;
    uvec doc;
    for (int i = 0; i < M; i++) {
        doc = as<uvec>(docs[i]);
        for (int j = 0; j < Ns(i); j++) {
            w = doc(j);
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    grad_THETA.row(i) += PHI_s(k, w-1) / ETA(i, w-1) * RHO(i,k) * 
                        ((double)(l==k) - RHO(i,l)) * (PSI.row(l) - THETA.row(i));
                }
            }
        }
    }

    // Gradients for PSI
    for (int i = 0; i < M; i++) {
        doc = as<uvec>(docs[i]);
        for (int j = 0; j < Ns(i); j++) {
            w = doc(j);
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    grad_PSI.row(k) += PHI_s(l, w-1) / ETA(i, w-1) * RHO(i,l) * 
                        ((double)(l==k) - RHO(i,k)) * (THETA.row(i) - PSI.row(k));
                }
            }
        }
    }

    // Gradient for PHI in simplex form.
    for (int i = 0; i < M; i++) {
        doc = as<uvec>(docs[i]);
        for (int j = 0; j < Ns(i); j++) {
            w = doc(j);
                for (int k = 0; k < K; k++) {
                grad_PHI(k, w-1) += RHO(i,k) / ETA(i, w-1);
            }
        }
    }

    // Propogate to real valued form if desired
    mat grad_GAMM(K, V-1, fill::zeros);
    for (int k = 0; k < K; k++) {
        for (int v = 0; v < V-1; v++) {
            for (int u = 0; u < V; u++) {
                grad_GAMM(k, v) += grad_PHI(k, u) * PHI_s(k,u) * 
                    (double(u==v) - PHI_s(k,v));
            }
        }
    }

    return List::create(Rcpp::Named("grad_THETA") = -grad_THETA,
                          Rcpp::Named("grad_PSI") = -grad_PSI,
                          Rcpp::Named("grad_PHI") = -grad_GAMM);
}
