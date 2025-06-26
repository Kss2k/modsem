#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "qml.h"
#include "mvnorm.h"

// [[Rcpp::export]]
arma::mat muQmlCpp(Rcpp::List m, int t, int ncores = 1) {
  ThreadSetter ts(ncores);                       // set threads

  int  numEta      = m["numEta"];
  int  numXi       = m["numXi"];
  arma::mat alpha  = m["alpha"];
  arma::mat beta0  = m["beta0"];
  arma::mat gammaXi= m["gammaXi"];
  arma::mat omegaXiXi = m["omegaXiXi"];
  arma::mat L1     = m["L1"];
  arma::mat L2     = m["L2"];
  arma::mat X      = m["x"];
  arma::mat U      = m["u"];
  arma::mat Sigma1 = m["Sigma1"];
  arma::mat Binv   = m["Binv"];
  arma::mat kronXi = m["kronXi"];

  arma::vec trOmegaSigma = traceOmegaSigma1(omegaXiXi * Sigma1, numEta);
  arma::mat Ey(t, numEta, arma::fill::none);

  const int lastColKOxx = numXi * numEta - 1,
            lastColBinv = numEta - 1;

  if (Binv.n_rows > static_cast<unsigned>(numEta)) {
    #ifdef _OPENMP
    #pragma omp parallel for if(ncores>1) schedule(static)
    #endif
    for (int i = 0; i < t; ++i) {
      int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;

      arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOxx);
      arma::mat Binv_t   = Binv  .submat(firstRow, 0, lastRow, lastColBinv);

      Ey.row(i) = ( Binv_t *
        ( trOmegaSigma + alpha
          + gammaXi * (beta0 + L1 * X.row(i).t())
          + kronXi_t * omegaXiXi * (beta0 + L1 * X.row(i).t()) )
        + L2 * U.row(i).t()
      ).t();
    }
  } else {
    #ifdef _OPENMP
    #pragma omp parallel for if(ncores>1) schedule(static)
    #endif
    for (int i = 0; i < t; ++i) {
      int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;

      arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOxx);

      Ey.row(i) = ( Binv *
        ( trOmegaSigma + alpha
          + gammaXi * (beta0 + L1 * X.row(i).t())
          + kronXi_t * omegaXiXi * (beta0 + L1 * X.row(i).t()) )
        + L2 * U.row(i).t()
      ).t();
    }
  }
  return Ey;
}


// [[Rcpp::export]]
arma::mat sigmaQmlCpp(Rcpp::List m, int t, int ncores = 1) {
  ThreadSetter ts(ncores);

  int  numEta           = m["numEta"];
  int  numXi            = m["numXi"];
  arma::mat gammaXi     = m["gammaXi"];
  arma::mat omegaXiXi   = m["omegaXiXi"];
  arma::mat L1          = m["L1"];
  arma::mat L2          = m["L2"];
  arma::mat X           = m["x"];
  arma::mat U           = m["u"];
  arma::mat Sigma1      = m["Sigma1"];
  arma::mat Sigma2Theta = m["Sigma2ThetaEpsilon"];
  arma::mat psi         = m["psi"];
  arma::mat Binv        = m["Binv"];
  arma::mat kronXi      = m["kronXi"];

  arma::mat varZ = varZCpp(omegaXiXi, Sigma1, numEta);
  arma::mat sigmaE(t * numEta, numEta, arma::fill::none);

  const int lastColSigma = numEta - 1,
            lastColKOxx  = numXi * numEta - 1;

  if (Binv.n_rows > static_cast<unsigned>(numEta)) {
    #ifdef _OPENMP
    #pragma omp parallel for if(ncores>1) schedule(static)
    #endif
    for (int i = 0; i < t; ++i) {
      int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;

      arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOxx);
      arma::mat Binv_t   = Binv  .submat(firstRow, 0, lastRow, lastColSigma);

      arma::mat Sigma2 = Binv_t * psi * Binv_t.t() + Sigma2Theta;
      sigmaE.submat(firstRow, 0, lastRow, lastColSigma) =
        ( Binv_t * (gammaXi + 2 * kronXi_t * omegaXiXi) ) * Sigma1 *
        ( Binv_t * (gammaXi + 2 * kronXi_t * omegaXiXi) ).t()
        + Sigma2 + Binv_t * varZ * Binv_t.t();
    }
  } else {                 // Binv is common to all time points
    varZ       = Binv * varZ * Binv.t();
    arma::mat Sigma2 = Binv * psi * Binv.t() + Sigma2Theta;
    #ifdef _OPENMP
    #pragma omp parallel for if(ncores>1) schedule(static)
    #endif
    for (int i = 0; i < t; ++i) {
      int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;

      arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOxx);

      sigmaE.submat(firstRow, 0, lastRow, lastColSigma) =
        ( Binv * (gammaXi + 2 * kronXi_t * omegaXiXi) ) * Sigma1 *
        ( Binv * (gammaXi + 2 * kronXi_t * omegaXiXi) ).t()
        + Sigma2 + varZ;
    }
  }
  return sigmaE;
}


// [[Rcpp::export]]
arma::mat calcKronXi(Rcpp::List m, int t, int ncores = 1) {
  ThreadSetter ts(ncores);

  int numEta  = m["numEta"];
  int numXi   = m["numXi"];
  arma::mat beta0 = m["beta0"];
  arma::mat L1    = m["L1"];
  arma::mat X     = m["x"];
  arma::mat Ie    = m["Ieta"];

  arma::mat out(t * numEta, numXi * numEta, arma::fill::none);

  #ifdef _OPENMP
  #pragma omp parallel for if(ncores>1) schedule(static)
  #endif
  for (int i = 0; i < t; ++i) {
    out.submat(i * numEta, 0,
               (i + 1) * numEta - 1, numXi * numEta - 1) =
      arma::kron(Ie, beta0.t() + X.row(i) * L1.t());
  }
  return out;
}


// [[Rcpp::export]]
arma::mat calcBinvCpp(Rcpp::List m, int t, int ncores = 1) {
  ThreadSetter ts(ncores);

  int numEta = m["numEta"],
      numXi  = m["numXi"],
      kOmega = m["kOmegaEta"];

  arma::mat gammaEta   = m["gammaEta"];
  arma::mat Ie         = m["Ieta"];
  arma::mat B0         = Ie - gammaEta;               // time-invariant part
  arma::mat omegaEtaXi = m["omegaEtaXi"];

  if (numEta == 1)           return Ie;
  if (kOmega == 0)           return arma::inv(B0);

  arma::mat kronXi = m["kronXi"];
  arma::mat Bout(t * numEta, numEta, arma::fill::none);

  const int lastColB   = numEta - 1,
            lastColKOx = numXi * numEta - 1;

  #ifdef _OPENMP
  #pragma omp parallel for if(ncores>1) schedule(static)
  #endif

  for (int i = 0; i < t; ++i) {
    int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;
    arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOx);
    Bout.submat(firstRow, 0, lastRow, lastColB) =
      arma::inv(B0 - kronXi_t * omegaEtaXi);
  }

  return Bout;
}


// [[Rcpp::export]]
arma::vec logNormalPdf(const arma::mat& X,
                       const arma::mat& mu,
                       const arma::vec& sigmaDiag,
                       int ncores = 1) {
  ThreadSetter ts(ncores);

  const double log2pi = std::log(2.0 * M_PI);
  arma::vec out(X.n_rows, arma::fill::none);

  #ifdef _OPENMP
  #pragma omp parallel for if(ncores>1) schedule(static)
  #endif
  for (arma::uword i = 0; i < X.n_rows; ++i) {
    double acc = 0.0;
    for (arma::uword j = 0; j < X.n_cols; ++j) {
      double diff = X(i, j) - mu(i, j);
      double sd   = sigmaDiag(j);
      acc += -0.5 * log2pi - std::log(sd) - 0.5 * (diff * diff) / (sd * sd);
    }
    out(i) = acc;
  }
  return out;
}


// [[Rcpp::export]]
arma::vec dnormCpp(const arma::vec& x, const arma::vec& mu, const arma::vec& sigma, 
                   int ncores = 1) {
  return logNormalPdf(x, mu, sigma, ncores);
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
