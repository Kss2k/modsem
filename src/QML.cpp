#include <RcppArmadillo.h>
#include "QML.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec muQmlCpp(Rcpp::List m, int t) {
  arma::mat alpha = Rcpp::as<arma::mat>(m["alpha"]);
  arma::mat gammaXi = Rcpp::as<arma::mat>(m["gammaXi"]);
  arma::mat omegaXiXi = Rcpp::as<arma::mat>(m["omegaXiXi"]);
  arma::mat l1 = Rcpp::as<arma::mat>(m["L1"]);
  arma::mat l2 = Rcpp::as<arma::mat>(m["L2"]);
  arma::mat x = Rcpp::as<arma::mat>(m["x"]);
  arma::mat u = Rcpp::as<arma::mat>(m["u"]);
  arma::mat Ey = arma::mat(t, 1);
  arma::mat Sigma1 = Rcpp::as<arma::mat>(m["Sigma1"]);

  double trOmegaSigma = arma::trace(omegaXiXi * Sigma1);
  for (int i = 0; i < t; i++) {
    Ey.row(i) = trOmegaSigma + alpha + gammaXi * l1 * x.row(i).t() + 
      x.row(i) * l1.t() * omegaXiXi * l1 * x.row(i).t() + l2 * u.row(i).t();
  }
  return Ey;
}


// [[Rcpp::export]]
arma::vec sigmaQmlCpp(Rcpp::List m, int t) {
  arma::mat gammaXi = Rcpp::as<arma::mat>(m["gammaXi"]);
  arma::mat omegaXiXi = Rcpp::as<arma::mat>(m["omegaXiXi"]);
  arma::mat l1 = Rcpp::as<arma::mat>(m["L1"]);
  arma::mat l2 = Rcpp::as<arma::mat>(m["L2"]);
  arma::mat x = Rcpp::as<arma::mat>(m["x"]);
  arma::mat u = Rcpp::as<arma::mat>(m["u"]);
  arma::mat Ey = arma::vec(t);
  arma::mat Sigma1 = Rcpp::as<arma::mat>(m["Sigma1"]);
  arma::mat Sigma2 = Rcpp::as<arma::mat>(m["Sigma2"]);
  
  arma::mat sigmaE = arma::mat(t, 1);
  double varZ = varZCpp(omegaXiXi, Sigma1);
  for (int i = 0; i < t; i++) {
    sigmaE.row(i) = (gammaXi + 2 * x.row(i) * l1.t() * omegaXiXi) * 
      Sigma1 * 
      (gammaXi + 2 * x.row(i) * l1.t() * omegaXiXi).t() + Sigma2  + varZ;
  }
  return sigmaE;
}


arma::vec logNormalPdf(const arma::vec& x, const arma::vec& mu, const arma::vec& sigma) {
    int n = x.n_elem;
    arma::vec result(n);

    // Constants
    double log_2pi = std::log(2.0 * M_PI);

    for (int i = 0; i < n; i++) {
        double diff = x(i) - mu(i);
        double sigma_sq = sigma(i) * sigma(i);

        // Log of the normal distribution PDF equation
        result(i) = -0.5 * log_2pi - std::log(sigma(i)) - 0.5 * (diff * diff) / sigma_sq;
    }

    return result;
}


// [[Rcpp::export]]
arma::vec dnormCpp(const arma::vec& x, const arma::vec& mu, const arma::vec& sigma) {
    return logNormalPdf(x, mu, sigma);
}


// [[Rcpp::export]]
double varZCpp(arma::mat Omega, arma::mat Sigma1) {

  int ds = Sigma1.n_rows;
  double varZ = 0;
  
  for (int i = 0; i < ds; i++) {
    for (int j = 0; j < ds; j++) {
      for (int k = 0; k < ds; k++) {
        for (int s = 0; s < ds; s++) {
          varZ += Omega(i, j) * Omega(k, s) * 
            (Sigma1(i, j) * Sigma1(k, s) + 
             Sigma1(i, k) * Sigma1(j, s) +
             Sigma1(i, s) * Sigma1(j, k));
        }
      }
    }
  }
  double trOmegaSigma1 = arma::trace(Omega * Sigma1);
  return varZ - trOmegaSigma1 * trOmegaSigma1;
}
