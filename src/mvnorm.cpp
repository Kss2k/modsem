// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "mvnorm.h"
#include <cfloat>


static double const log2pi = std::log(2.0 * M_PI);


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}


// [[Rcpp::export]]
arma::vec rep_dmvnorm(const arma::mat &x, 
                      const arma::mat &expected, 
                      const arma::mat &sigma, const int t) {
  arma::vec out = arma::vec(t);
  int ncolSigma = sigma.n_cols;
  int firstRow, lastRow, firstCol, lastCol;
  for (int i = 0; i < t; i++) {
    firstRow = i * ncolSigma;
    lastRow = (i + 1) * ncolSigma - 1;
    firstCol = 0;
    lastCol = ncolSigma - 1;

    out(i) = dmvnrm_arma_mc(x.row(i), expected.row(i), 
        sigma.submat(firstRow, firstCol, lastRow, lastCol), true)(0);
  }
  return out;
}


// [[Rcpp::export]]
arma::vec dmvnrm_arma_mc(arma::mat const &x,  
                         arma::rowvec const &mean,  
                         arma::mat const &sigma, 
                         bool const logd = true) {
    using arma::uword;
    uword const n = x.n_rows, 
             xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
                constants = -(double)xdim/2.0 * log2pi, 
              other_terms = rootisum + constants;
    
    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);   
        out(i) = other_terms - 0.5 * arma::dot(z, z);     
    }  
      
    if (logd)
      return out;
    return exp(out);
}


// [[Rcpp::export]]
double totalDmvnWeightedCpp(const arma::vec& mu,
                            const arma::mat& sigma,
                            const arma::vec& nu,
                            const arma::mat& S,
                            double tgamma,
                            int n,
                            int d) {
  if (!sigma.is_finite()) return NA_REAL;

  if (arma::any(arma::diagvec(sigma) <= 0)) return NA_REAL;

  arma::mat L;
  if (!arma::chol(L, sigma, "lower")) return NA_REAL;

  double log_det_sigma = 2.0 * arma::sum(log(L.diag()));
  arma::mat sigma_inv = arma::inv_sympd(sigma);
  if (!sigma_inv.is_finite()) return NA_REAL;

  arma::vec diff = nu - mu;
  double mahalanobis_term = tgamma * arma::as_scalar(diff.t() * sigma_inv * diff);

  double trace_term = arma::accu(sigma_inv % S);

  double log_likelihood = -0.5 * (
    tgamma * d * std::log(2.0 * M_PI) +
    tgamma * log_det_sigma +
    trace_term +
    mahalanobis_term
  );

  return log_likelihood;
}


arma::vec dmvnorm_log(const arma::mat& X,
                             const arma::vec& mu,
                             const arma::mat& L) {
    arma::mat Z = (X.each_row() - mu.t()) *
                  arma::inv(trimatu(L.t()));

    arma::vec quad = arma::sum(arma::square(Z), 1); 
    double logdet  = 2.0 * arma::sum(log(L.diag()));

    double cst = -0.5 * ( X.n_cols * std::log(2.0 * M_PI) + logdet );
    return cst - 0.5 * quad;                           
}
