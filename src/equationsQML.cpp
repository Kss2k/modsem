#include <RcppArmadillo.h>
#include "QML.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat muQmlCpp(Rcpp::List m, int t) {
  int numEta = Rcpp::as<int>(m["numEta"]);
  int numXi = Rcpp::as<int>(m["numXi"]);
  arma::mat alpha = Rcpp::as<arma::mat>(m["alpha"]);
  arma::mat beta0 = Rcpp::as<arma::mat>(m["beta0"]);
  arma::mat gammaXi = Rcpp::as<arma::mat>(m["gammaXi"]);
  arma::mat omegaXiXi = Rcpp::as<arma::mat>(m["omegaXiXi"]);
  arma::mat l1 = Rcpp::as<arma::mat>(m["L1"]); // L1 refers to L1 cache
  arma::mat l2 = Rcpp::as<arma::mat>(m["L2"]); // L2 refers to L2 cache
  arma::mat x = Rcpp::as<arma::mat>(m["x"]);
  arma::mat u = Rcpp::as<arma::mat>(m["u"]);
  arma::mat Ey = arma::mat(t, numEta);
  arma::mat Sigma1 = Rcpp::as<arma::mat>(m["Sigma1"]);
  arma::mat Ie = Rcpp::as<arma::mat>(m["Ieta"]);
  arma::mat Binv = Rcpp::as<arma::mat>(m["Binv"]);
  arma::vec trOmegaSigma = traceOmegaSigma1(omegaXiXi * Sigma1, numEta);
  arma::mat kronXi = Rcpp::as<arma::mat>(m["kronXi"]); 

  arma::mat kronXi_t;
  
  int firstRow, lastRow, firstCol = 0,
      lastColKOxx = numXi * numEta - 1, lastColBinv = numEta - 1;

  if (int(Binv.n_rows) > numEta) {
    arma::mat Binv_t;
    for (int i = 0; i < t; i++) {
      firstRow = i * numEta;
      lastRow = (i + 1) * numEta - 1;

      kronXi_t = 
        kronXi.submat(firstRow, firstCol, lastRow, lastColKOxx);
      Binv_t = Binv.submat(firstRow, firstCol, lastRow, lastColBinv);

      Ey.row(i) = (Binv_t * (trOmegaSigma + alpha + gammaXi * (beta0 + l1 * x.row(i).t()) + 
        kronXi_t * omegaXiXi * (beta0 + l1 * x.row(i).t())) + l2 * u.row(i).t()).t();
    }
  } else {
    for (int i = 0; i < t; i++) {
      firstRow = i * numEta;
      lastRow = (i + 1) * numEta - 1;

      kronXi_t = 
        kronXi.submat(firstRow, firstCol, lastRow, lastColKOxx);

      Ey.row(i) = (Binv * (trOmegaSigma + alpha + gammaXi * (beta0 + l1 * x.row(i).t()) + 
            kronXi_t * omegaXiXi * (beta0 + l1 * x.row(i).t())) + l2 * u.row(i).t()).t();
    }
  }
  return Ey;
}


// [[Rcpp::export]]
arma::mat sigmaQmlCpp(Rcpp::List m, int t) {
  int numEta = Rcpp::as<int>(m["numEta"]);
  int numXi = Rcpp::as<int>(m["numXi"]);
  arma::mat gammaXi = Rcpp::as<arma::mat>(m["gammaXi"]);
  arma::mat omegaXiXi = Rcpp::as<arma::mat>(m["omegaXiXi"]);
  arma::mat l1 = Rcpp::as<arma::mat>(m["L1"]);
  arma::mat l2 = Rcpp::as<arma::mat>(m["L2"]);
  arma::mat x = Rcpp::as<arma::mat>(m["x"]);
  arma::mat u = Rcpp::as<arma::mat>(m["u"]);
  arma::mat Sigma1 = Rcpp::as<arma::mat>(m["Sigma1"]);
  arma::mat Sigma2ThetaEpsilon = Rcpp::as<arma::mat>(m["Sigma2ThetaEpsilon"]);
  arma::mat psi = Rcpp::as<arma::mat>(m["psi"]);
  arma::mat Ie = Rcpp::as<arma::mat>(m["Ieta"]);
  arma::mat sigmaE = arma::mat(t * numEta, numEta);
  arma::mat Binv = Rcpp::as<arma::mat>(m["Binv"]);
  arma::mat varZ = varZCpp(omegaXiXi, Sigma1, numEta); 
  arma::mat  kronXi = Rcpp::as<arma::mat>(m["kronXi"]); 
  
  int firstRow, lastRow, firstCol = 0, lastColSigmaE = numEta - 1, 
      lastColKOxx = numXi * numEta - 1;
  arma::mat kronXi_t;
  arma::mat Sigma2;
  if (int(Binv.n_rows) > numEta) {
    arma::mat Binv_t;
    for (int i = 0; i < t; i++) {
      firstRow = i * numEta;
      lastRow = (i + 1) * numEta - 1;
      
      kronXi_t = 
        kronXi.submat(firstRow, firstCol, lastRow, lastColKOxx);
      Binv_t = Binv.submat(firstRow, firstCol, lastRow, lastColSigmaE);

      Sigma2 = Binv_t * psi * Binv_t.t() + Sigma2ThetaEpsilon;
      sigmaE.submat(firstRow, firstCol, lastRow, lastColSigmaE) = 
        (Binv_t * (gammaXi + 2 * kronXi_t * omegaXiXi)) * Sigma1 * 
        (Binv_t * (gammaXi + 2 * kronXi_t * omegaXiXi)).t() + Sigma2 + 
        Binv_t * varZ * Binv_t.t();
    }
  } else {
    varZ = Binv * varZ * Binv.t();
    Sigma2 = Binv * psi * Binv.t() + Sigma2ThetaEpsilon;
    for (int i = 0; i < t; i++) {
      firstRow = i * numEta;
      lastRow = (i + 1) * numEta - 1;

      kronXi_t = 
        kronXi.submat(firstRow, firstCol, lastRow, lastColKOxx);

      sigmaE.submat(firstRow, firstCol, lastRow, lastColSigmaE) = 
        (Binv * (gammaXi + 2 * kronXi_t * omegaXiXi)) * Sigma1 * 
        (Binv * (gammaXi + 2 * kronXi_t * omegaXiXi)).t() + Sigma2 + varZ;
    }
  }

  return sigmaE;
}


// [[Rcpp::export]]
arma::mat calcKronXi(Rcpp::List m, int t) {
  int numEta = Rcpp::as<int>(m["numEta"]);  
  int numXi = Rcpp::as<int>(m["numXi"]);
  arma::mat beta0 = Rcpp::as<arma::mat>(m["beta0"]);
  arma::mat omegaXiXi = Rcpp::as<arma::mat>(m["omegaXiXi"]);
  arma::mat l1 = Rcpp::as<arma::mat>(m["L1"]); 
  arma::mat x = Rcpp::as<arma::mat>(m["x"]);
  arma::mat Ie = Rcpp::as<arma::mat>(m["Ieta"]);
  // dimensions of a single kronecker product (A (kron) B) = (m x n) (kron) (p x q) = (mp x nq)
  // in this case: (numEta x numEta) (kron) (1 x numXi) = (numEta x numEta * numXi) 
  arma::mat out = arma::mat(t * numEta, numXi * numEta);

  for (int i = 0; i < t; i++) {
    out.submat(i * numEta, 0, (i + 1) * numEta - 1, numXi * numEta - 1) = 
      arma::kron(Ie, beta0.t() + x.row(i) * l1.t());
  }
  return out;
}


// [[Rcpp::export]]
arma::mat calcBinvCpp(Rcpp::List m, int t) {
  int numEta = Rcpp::as<int>(m["numEta"]); 
  int numXi = Rcpp::as<int>(m["numXi"]);
  int kOmegaEta = Rcpp::as<int>(m["kOmegaEta"]);
  arma::mat gammaEta = Rcpp::as<arma::mat>(m["gammaEta"]);

  arma::mat Ie = Rcpp::as<arma::mat>(m["Ieta"]); 
  arma::mat B = Ie - gammaEta;
  arma::mat omegaEtaXi = Rcpp::as<arma::mat>(m["omegaEtaXi"]);

  if (numEta == 1) return Ie;
  else if (kOmegaEta == 0) return arma::inv(B);

  arma::mat kronXi = Rcpp::as<arma::mat>(m["kronXi"]);
  arma::mat B_t = arma::mat(t * numEta, numEta);

  int firstRow, lastRow, firstCol = 0, lastColB = numEta - 1, 
      lastColKOxx = numXi * numEta - 1;
  arma::mat kronXi_t;

  for (int i = 0; i < t; i++) {
    firstRow = i * numEta;
    lastRow = (i + 1) * numEta - 1;
    
    kronXi_t = 
      kronXi.submat(firstRow, firstCol, lastRow, lastColKOxx);

    B_t.submat(firstRow, firstCol, lastRow, lastColB) = 
      arma::inv(B - kronXi_t * omegaEtaXi);
  }

  return B_t;
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
