#include <RcppArmadillo.h>
#include "lms.h"
#include "mvnorm.h"
#include <float.h>

// [[Rcpp::depends(RcppArmadillo)]]


// Deprecated will remove soon...
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


// Deprecated will remove soon...
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


inline arma::mat make_Oi(unsigned k, unsigned numXis) {
  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));
  return Oi;
}


inline arma::vec make_zvec(unsigned k, unsigned numXis, const arma::vec& z) {
  if (k > 0) return arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else return arma::zeros<arma::vec>(numXis);
}


struct LMSModel {
  arma::mat A, Oxx, Oex, Ie, lY, lX, tY, tX, Gx, Ge,
    a, beta0, Psi, d, e;
  unsigned  k       = 0;
  unsigned  numXis  = 0;

  explicit LMSModel(const Rcpp::List& modFilled) {

    Rcpp::List matrices = modFilled["matrices"];
    Rcpp::List info     = modFilled["info"];
    Rcpp::List quad     = modFilled["quad"];

    k       = Rcpp::as<unsigned>(quad["k"]);
    numXis  = Rcpp::as<unsigned>(info["numXis"]);

    // one-liners, no loops
    A      = Rcpp::as<arma::mat>(matrices["A"]);
    Oxx    = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
    Oex    = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
    Ie     = Rcpp::as<arma::mat>(matrices["Ieta"]);
    lY     = Rcpp::as<arma::mat>(matrices["lambdaY"]);
    lX     = Rcpp::as<arma::mat>(matrices["lambdaX"]);
    tY     = Rcpp::as<arma::mat>(matrices["tauY"]);
    tX     = Rcpp::as<arma::mat>(matrices["tauX"]);
    Gx     = Rcpp::as<arma::mat>(matrices["gammaXi"]);
    Ge     = Rcpp::as<arma::mat>(matrices["gammaEta"]);
    a      = Rcpp::as<arma::mat>(matrices["alpha"]);
    beta0  = Rcpp::as<arma::mat>(matrices["beta0"]);
    Psi    = Rcpp::as<arma::mat>(matrices["psi"]);
    d      = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
    e      = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);
  }

  arma::vec mu(const arma::vec& z) const {
    arma::vec zVec = make_zvec(k, numXis, z);
    arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

    arma::mat Binv;
    if (Ie.n_cols == 1) Binv = arma::mat(Ie);
    else Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    arma::vec muX = tX + lX * (beta0 + A * zVec);
    arma::vec muY = tY + 
      lY * (Binv * (a + 
            Gx * (beta0 + A * zVec) + 
            kronZ.t() * Oxx * (beta0 + A * zVec)));
    return arma::join_cols(muX, muY);
  }


  arma::mat Sigma(const arma::vec& z) const {
    arma::vec zVec = make_zvec(k, numXis, z);

    arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

    arma::mat Binv;
    if (Ie.n_cols == 1) Binv = arma::mat(Ie);
    else Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    arma::mat Oi = make_Oi(k, numXis);
    arma::mat Sxx = lX * A * Oi * A.t() * lX.t() + d;
    arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
    arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
    arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() + 
      lY * (Binv * Psi * Binv.t()) * lY.t() + e;

    return arma::join_cols(
        arma::join_rows(Sxx, Sxy), 
        arma::join_rows(Sxy.t(), Syy)
        );
  }
};


// [[Rcpp::export]]
double completeLogLikLmsCpp(Rcpp::List modelR, Rcpp::List P, Rcpp::List quad) {
  const LMSModel model = LMSModel(modelR);
  const arma::mat Z = Rcpp::as<arma::mat>(quad["n"]).t(); // transpose so we can use column order vectors

  Rcpp::List info   = modelR["info"];
  Rcpp::List TGamma = P["tgamma"];
  Rcpp::List Means  = P["mean"];
  Rcpp::List Covs   = P["cov"];

  const int n = Rcpp::as<int>(info["N"]);
  const int d = Rcpp::as<int>(info["ncol"]);

  double ll = 0;
  for (std::size_t i = 0; i < Z.n_cols; i++) {
    const double tgamma = Rcpp::as<double>(TGamma[i]);
    if (tgamma <= DBL_MIN) continue;

    const arma::vec z = Z.col(i);
    const arma::vec nu = Rcpp::as<arma::vec>(Means[i]);
    const arma::mat S = Rcpp::as<arma::mat>(Covs[i]);
    const arma::mat mu = model.mu(z);
    const arma::mat Sigma = model.Sigma(z);

    ll += totalDmvnWeightedCpp(mu, Sigma, nu, S, tgamma, n, d);
  }

  return ll;
}

inline double
logLikFromModel(const LMSModel&            M,
                const arma::mat&           V,          // nodes × k   (quad$n)
                const Rcpp::List&          Mean,       // P$mean
                const Rcpp::List&          Cov,        // P$cov
                const Rcpp::List&          TGamma,     // P$tgamma
                const int                  n,          // # rows data
                const int                  d)          // # manifests
{
    const std::size_t J = V.n_rows;
    double ll = 0.0;

    for (std::size_t j = 0; j < J; ++j) {

        const double tg = Rcpp::as<double>(TGamma[j]);
        if (tg <= DBL_MIN) continue;

        const arma::vec z   = V.row(j).t();
        const arma::vec nu  = Rcpp::as<arma::vec>(Mean[j]);
        const arma::mat S   = Rcpp::as<arma::mat>(Cov [j]);

        const arma::vec mu  = M.mu   (z);
        const arma::mat Sig = M.Sigma(z);

        ll += totalDmvnWeightedCpp(mu, Sig, nu, S, tg, n, d);
    }
    return ll;
}


// [[Rcpp::export]]
arma::vec gradLogLikLmsCpp(const Rcpp::List& modelR,      // *filled* model
                           const Rcpp::List& P,           // E-step cache
                           const arma::uvec& block,       // free-param blocks
                           const arma::uvec& row,
                           const arma::uvec& col,
                           double            eps = 1e-6) { 
  const LMSModel Mbase(modelR);           // build C++ struct once

  const arma::mat V       = Rcpp::as<arma::mat>(P["V"]);
  const Rcpp::List Mean   = P["mean"];
  const Rcpp::List Cov    = P["cov"];
  const Rcpp::List TGamma = P["tgamma"];

  const Rcpp::List info   = modelR["info"];
  const int n             = Rcpp::as<int>(info["N"]);
  const int d             = Rcpp::as<int>(info["ncol"]);

  // base likelihood
  double f0 = logLikFromModel(Mbase, V, Mean, Cov, TGamma, n, d);

  // will mutate *one* entry of these matrices in place
  LMSModel M = Mbase;                     // lightweight copy (pointers)

  const std::size_t nPar = block.n_elem;
  arma::vec grad(nPar);

  for (std::size_t k = 0; k < nPar; ++k) {

    double* target = nullptr;

    switch (block[k]) {
      case 0:  target = &M.lX   (row[k], col[k]); break;
      case 1:  target = &M.lY   (row[k], col[k]); break;
      case 2:  target = &M.tX   (row[k], col[k]); break;
      case 3:  target = &M.tY   (row[k], col[k]); break;
      case 4:  target = &M.d    (row[k], col[k]); break;
      case 5:  target = &M.e    (row[k], col[k]); break;
      case 6:  target = &M.A    (row[k], col[k]); break;
      case 7:  target = &M.Psi  (row[k], col[k]); break;
      case 8:  target = &M.a    (row[k], col[k]); break;
      case 9:  target = &M.beta0(row[k], col[k]); break;
      case 10: target = &M.Gx   (row[k], col[k]); break;
      case 11: target = &M.Ge   (row[k], col[k]); break;
      case 12: target = &M.Oxx  (row[k], col[k]); break;
      case 13: target = &M.Oex  (row[k], col[k]); break;
      default: Rcpp::stop("unknown block");
    }

    const double old = *target;
    *target += eps;

    const double f1 = logLikFromModel(M, V, Mean, Cov, TGamma, n, d);

    grad[k] = (f1 - f0) / eps;

    *target = old;                         // restore for next param
  }

  return grad;
}
