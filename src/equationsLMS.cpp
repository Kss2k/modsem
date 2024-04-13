#include <RcppArmadillo.h>
#include "LMS.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec muLmsCpp(Rcpp::List model, arma::vec z) {
  Rcpp::List matrices = model["matrices"];
  Rcpp::List info = model["info"];
  Rcpp::List quad = model["quad"];
  int numEtas = Rcpp::as<int>(info["numEtas"]);
  int numXis = Rcpp::as<int>(info["numXis"]);
  int k = Rcpp::as<int>(quad["k"]);
  arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  arma::mat subA = A.submat(0, 0, numXis - 1, numXis - 1);
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
  
  arma::vec collapseEta = arma::ones<arma::vec>(numEtas);
  arma::mat cOex = Rcpp::as<arma::mat>(matrices["selectionMatrixOmegaEtaXi"]);
  arma::vec zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  arma::mat zMat = zToMatrix(zVec, numEtas);
  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - (zMat.t() * A.t() * Oex) * cOex);
  }
  arma::vec muX = tX + lX * subA * zVec;
  arma::vec muY = tY + 
    lY * (Binv * (a + 
            Gx * subA * zVec + 
            (zMat.t() * A.t() * Oxx * A * zMat) * 
            collapseEta));
  return arma::join_cols(muX, muY);
}


// [[Rcpp::export]]
arma::mat sigmaLmsCpp(Rcpp::List model, arma::vec z) {
  Rcpp::List matrices = model["matrices"];
  Rcpp::List info = model["info"];
  Rcpp::List quad = model["quad"];
  int numEtas = Rcpp::as<int>(info["numEtas"]);
  int numXis = Rcpp::as<int>(info["numXis"]);
  int k = Rcpp::as<int>(quad["k"]);
  arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  arma::mat subA = A.submat(0, 0, numXis - 1, numXis - 1);
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
  arma::mat Psi = Rcpp::as<arma::mat>(matrices["psi"]); 
  arma::mat d = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  arma::mat e = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);
  arma::vec collapseEta = arma::ones<arma::vec>(numEtas);
  arma::mat cOex = Rcpp::as<arma::mat>(matrices["selectionMatrixOmegaEtaXi"]);
  arma::mat cOxx = Rcpp::as<arma::mat>(matrices["selectionMatrixOmega"]);
  arma::vec zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  arma::mat zMat = zToMatrix(zVec, numEtas);
  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - (zMat.t() * A.t() * Oex) * cOex);
  }
  
  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));
  arma::mat Sxx = lX * subA * Oi * subA.t() * lX.t() + d;
  arma::mat Eta = Binv * (Gx * subA + (zMat.t() * A.t() * Oxx * A) * cOxx);
  arma::mat Sxy = lX * (subA * Oi * Eta.t()) * lY.t();
  arma::mat Syy = lY * (Eta * Oi * Eta.t()) * lY.t() + 
    lY * (Binv * Psi * Binv.t()) * lY.t() + e;
  return arma::join_cols(arma::join_rows(Sxx, Sxy), arma::join_rows(Sxy.t(), Syy));
}


// [[Rcpp::export]]
arma::mat zToMatrix(arma::vec z, int numEtas) {
  arma::mat mat = arma::zeros<arma::mat>(numEtas * z.n_elem, numEtas);
  for (int i = 0; i < numEtas; i++) {
    mat.submat(i * z.n_elem, i, (i + 1) * z.n_elem - 1, i) = z;
  }
  return mat;
}
