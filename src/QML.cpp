#include <RcppArmadillo.h>
#include "QML.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat muQmlCpp(Rcpp::List m, int t) {
  int numEta = Rcpp::as<int>(m["numEta"]);
  arma::mat alpha = Rcpp::as<arma::mat>(m["alpha"]);
  arma::mat gammaXi = Rcpp::as<arma::mat>(m["gammaXi"]);
  arma::mat omegaXiXi = Rcpp::as<arma::mat>(m["omegaXiXi"]);
  arma::mat l1 = Rcpp::as<arma::mat>(m["L1"]); // L1 refers to L1 cache
  arma::mat l2 = Rcpp::as<arma::mat>(m["L2"]); // L2 refers to L2 cache
  arma::mat x = Rcpp::as<arma::mat>(m["x"]);
  arma::mat u = Rcpp::as<arma::mat>(m["u"]);
  arma::mat Ey = arma::mat(t, numEta);
  arma::mat Sigma1 = Rcpp::as<arma::mat>(m["Sigma1"]);
  arma::mat Ie = Rcpp::as<arma::mat>(m["Ieta"]);
  arma::vec trOmegaSigma = traceOmegaSigma1(omegaXiXi * Sigma1, numEta);

  for (int i = 0; i < t; i++) {
    Ey.row(i) = (trOmegaSigma + alpha + gammaXi * l1 * x.row(i).t() + 
      arma::kron(Ie, x.row(i) * l1.t()) * omegaXiXi * l1 * x.row(i).t() + l2 * u.row(i).t()).t();
  }
  return Ey;
}


// [[Rcpp::export]]
arma::mat sigmaQmlCpp(Rcpp::List m, int t) {
  int numEta = Rcpp::as<int>(m["numEta"]);
  arma::mat gammaXi = Rcpp::as<arma::mat>(m["gammaXi"]);
  arma::mat omegaXiXi = Rcpp::as<arma::mat>(m["omegaXiXi"]);
  arma::mat l1 = Rcpp::as<arma::mat>(m["L1"]);
  arma::mat l2 = Rcpp::as<arma::mat>(m["L2"]);
  arma::mat x = Rcpp::as<arma::mat>(m["x"]);
  arma::mat u = Rcpp::as<arma::mat>(m["u"]);
  arma::mat Ey = arma::vec(t);
  arma::mat Sigma1 = Rcpp::as<arma::mat>(m["Sigma1"]);
  arma::mat Sigma2 = Rcpp::as<arma::mat>(m["Sigma2"]);
  arma::mat kronXiOmega; 
  arma::mat Ie = Rcpp::as<arma::mat>(m["Ieta"]);
  arma::mat sumVec = Rcpp::as<arma::vec>(m["sumVec"]);
  arma::mat sigmaE = arma::mat(t, numEta);
  arma::mat varZ = varZCpp(omegaXiXi, Sigma1, numEta); // uneccessary to call this twice?
  
  for (int i = 0; i < t; i++) {
    kronXiOmega = arma::kron(Ie, x.row(i) * l1.t()) * omegaXiXi; 
    sigmaE.row(i) = 
      (arma::diagmat((gammaXi + 2 * kronXiOmega) * Sigma1 * 
      (gammaXi + 2 * kronXiOmega).t() + Sigma2  + varZ) * sumVec).t();
  }
  return sigmaE;
}


arma::vec logNormalPdf(const arma::vec& x, const arma::vec& mu, const arma::mat& sigma) {

  int n = x.n_elem;
  arma::vec result = arma::zeros<arma::vec>(n);

  // Constants
  double log_2pi = std::log(2.0 * M_PI);
  for (int j = 0; j < int(sigma.n_cols); j++) {
    for (int i = 0; i < n; i++) {
      double diff = x(i, j) - mu(i, j);
      double sigma_sq = sigma(i, j) * sigma(i, j);

      // Log of the normal distribution PDF equation
      result(i) += -0.5 * log_2pi - std::log(sigma(i, j)) - 0.5 * (diff * diff) / sigma_sq;
    }
  }
  return result;
}


// [[Rcpp::export]]
arma::vec dnormCpp(const arma::vec& x, const arma::vec& mu, const arma::vec& sigma) {
    return logNormalPdf(x, mu, sigma);
}


// [[Rcpp::export]]
arma::mat varZCpp(arma::mat Omega, arma::mat Sigma1, int numEta) {
  arma::mat varZ = arma::mat(numEta, numEta);
  int subRows = Omega.n_rows / numEta; 
  for (int i = 0; i < numEta; i++) {
    varZ(i, i) = varZSubOmega(Omega.submat(i * subRows, 0,
          (i + 1) * subRows - 1, (Omega.n_cols - 1)), Sigma1);
  }
  return varZ;
}


double varZSubOmega(arma::mat Omega, arma::mat Sigma1) {

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


arma::vec traceOmegaSigma1(const arma::mat OmegaSigma1, const int numEta) {
  arma::vec trace = arma::vec(numEta);
  int subRows = OmegaSigma1.n_rows / numEta;
  for (int i = 0; i < numEta; i++) {
    for (int j = 0; j < int(OmegaSigma1.n_cols); j++) {
      trace(i) += OmegaSigma1(i * subRows + j, j);
    } 
  }
  return trace;
}
