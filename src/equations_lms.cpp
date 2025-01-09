#include <RcppArmadillo.h>
#include "lms.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec muLmsCpp(Rcpp::List model, arma::vec z) {
  Rcpp::List matrices = model["matrices"];
  Rcpp::List info = model["info"];
  Rcpp::List quad = model["quad"];
  int numXis = Rcpp::as<int>(info["numXis"]);
  int k = Rcpp::as<int>(quad["k"]);
  arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  arma::mat Ie = Rcpp::as<arma::mat>(matrices["Ieta"]);
  arma::mat lY = Rcpp::as<arma::mat>(matrices["lambdaY"]);
  arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  arma::mat tY = Rcpp::as<arma::mat>(matrices["tauY"]);
  arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);
  
  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else zVec = arma::zeros<arma::vec>(numXis);
  arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);
  }
  arma::vec muX = tX + lX * (beta0 + A * zVec);
  arma::vec muY = tY + 
    lY * (Binv * (a + 
          Gx * (beta0 + A * zVec) + 
          kronZ.t() * Oxx * (beta0 + A * zVec)));
  return arma::join_cols(muX, muY);
}


// [[Rcpp::export]]
arma::mat sigmaLmsCpp(Rcpp::List model, arma::vec z) {
  Rcpp::List matrices = model["matrices"];
  Rcpp::List info = model["info"];
  Rcpp::List quad = model["quad"];
  int numXis = Rcpp::as<int>(info["numXis"]);
  int k = Rcpp::as<int>(quad["k"]);
  arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  arma::mat Ie = Rcpp::as<arma::mat>(matrices["Ieta"]);
  arma::mat lY = Rcpp::as<arma::mat>(matrices["lambdaY"]);
  arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  arma::mat tY = Rcpp::as<arma::mat>(matrices["tauY"]);
  arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);
  arma::mat Psi = Rcpp::as<arma::mat>(matrices["psi"]); 
  arma::mat d = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  arma::mat e = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else zVec = arma::zeros<arma::vec>(numXis);
  arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);
  }
  
  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));
  arma::mat Sxx = lX * A * Oi * A.t() * lX.t() + d;
  arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
  arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
  arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() + 
    lY * (Binv * Psi * Binv.t()) * lY.t() + e;
  return arma::join_cols(arma::join_rows(Sxx, Sxy), arma::join_rows(Sxy.t(), Syy));
}
