/** src/exp_nlpost.cpp Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.08.2018  */

#include <Rcpp.h>
using namespace Rcpp;

// Note: parameterized in terms of phi = precision = 1/variance = 1/sqrt(std)
double logdnorm(double x, double mu, double phi) {
    // TODO: Remove normalizing constant.
    return 0.5 * log(phi / (2.0 * M_PI)) - 0.5 * pow(x - mu, 2) * phi;
}

// Numerically stable softmax function.
// [[Rcpp::export]]
NumericVector softmaxC(NumericVector x) {
    return exp(x - max(x) - log(sum(exp(x - max(x)))));

}

// [[Rcpp::export]]
double exp_nlpostC(List Z_exp, NumericMatrix PHI, NumericMatrix THETA, 
        NumericMatrix PSI, List docs, NumericVector Ns, double eta, double gamma, 
        double beta) {

    int K = PHI.nrow();
    int M = THETA.nrow();
    int P = THETA.ncol();

    // -- Prior Distribution -- 
    double ll = 0;
    // Topic by Word matrix (dirichlet prior)
    for (int k = 0; k < K; k++) {
        ll += (eta - 1) * sum(log(PHI(k,_)));
    }

    // Topic Locations
    for (int k = 0; k < K; k++) {
        for (int p = 0; p < P; p++) {
            ll += logdnorm(PSI(k, p), 0, beta);
        }
    }
    // Doc Locations
    for (int m = 0; m < M; m++) {
        for (int p = 0; p < P; p++) {
            ll += logdnorm(THETA(m, p), 0, gamma);
        }
    }

    // -- Expected Likelihood -- 
    // Get probability of topic in each doc
    NumericMatrix RHO(M, K);
    NumericVector dists(K);
    for (int m = 0; m < M; m++) {
        for (int k = 0; k < K; k++) {
            dists(k) = sum(pow(THETA(m,_) - PSI(k,_), 2));
        }
        RHO(m,_) = softmaxC(-0.5 * dists);
    }

    // Sum up probs for each doc
    for (int m = 0; m < M; m++) {
        IntegerVector doc = as<IntegerVector>(docs[m]);
        NumericMatrix Z = as<NumericMatrix>(Z_exp[m]);

        for (int n = 0; n < Ns(m); n++) {
            for (int k = 0; k < K; k++) {
                ll += Z(n,k) * log(RHO(m, k) * PHI(k, doc(n) - 1));
            }
        }
    }

    return(-ll);
}
