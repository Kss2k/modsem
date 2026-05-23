#include <RcppArmadillo.h>
#include "mvnorm.h"

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat diagPartitioned(const arma::mat X, const arma::mat Y) {
  const int nx = X.n_rows;
  const int mx = X.n_cols;
  const int ny = Y.n_rows;
  const int my = Y.n_cols;

  if      (ny <= 0 || my <= 0) return(X);
  else if (nx <= 0 || mx <= 0) return(Y);

  arma::mat Z(nx + ny, mx + my, arma::fill::zeros);
  Z.submat(0, 0, nx - 1, mx - 1) = X;
  Z.submat(nx, mx, nx + ny - 1, mx + my - 1) = Y;

  return Z;
}

struct ModelMatrices {
  // Measurement Model
  arma::mat Lambda;
  arma::mat Theta;
  arma::vec tau;

  // Structural Model
  arma::mat Psi;
  arma::mat GammaXi;
  arma::mat GammaEta;
  arma::mat OmegaXiXi;
  arma::mat OmegaEtaXi;
  arma::vec beta0;
  arma::vec alpha;
  arma::mat Ie;

  ModelMatrices(const Rcpp::List matrices) {
    // Measurement Model
    const arma::mat LambdaX      = Rcpp::as<arma::mat>(matrices["lambdaX"]);
    const arma::mat LambdaY      = Rcpp::as<arma::mat>(matrices["lambdaY"]);
    const arma::mat ThetaDelta   = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
    const arma::mat ThetaEpsilon = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);
    const arma::vec tauX         = Rcpp::as<arma::vec>(matrices["tauX"]);
    const arma::vec tauY         = Rcpp::as<arma::vec>(matrices["tauY"]);

    Lambda = diagPartitioned(LambdaX, LambdaY);
    Theta  = diagPartitioned(ThetaDelta, ThetaEpsilon);
    tau    = arma::join_cols(tauX, tauY);

    // Structural Model
    const arma::mat PhiX = Rcpp::as<arma::mat>(matrices["phi"]);
    const arma::mat PsiY = Rcpp::as<arma::mat>(matrices["psi"]);

    Psi        = diagPartitioned(PhiX, PsiY);
    GammaXi    = Rcpp::as<arma::mat>(matrices["gammaXi"]);
    GammaEta   = Rcpp::as<arma::mat>(matrices["gammaEta"]);
    OmegaXiXi  = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
    OmegaEtaXi = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
    beta0      = Rcpp::as<arma::vec>(matrices["beta0"]);
    alpha      = Rcpp::as<arma::vec>(matrices["alpha"]);
    Ie         = Rcpp::as<arma::mat>(matrices["Ieta"]);
  }
};


// [[Rcpp::export]]
SEXP modelMatrixCacheCpp(Rcpp::List matrices) {
  ModelMatrices* M = new ModelMatrices(matrices);
  return Rcpp::XPtr<ModelMatrices>(M, true);
}


// [[Rcpp::export]]
arma::vec impliedEtaFromZetaCpp(
    const arma::vec zeta,
    const SEXP xptr
) {
  const Rcpp::XPtr<ModelMatrices> M(xptr);

  // dim
  const int k = M->GammaEta.n_cols;
  const int m = M->GammaXi.n_cols;

  if (k + m != (int)zeta.n_elem)
    Rcpp::stop("Wrong dimension of zeta!");

  // Split Zeta into exogenous and endogenous blocks
  const arma::mat xi = M->beta0 + zeta.subvec(0, m - 1);
  const arma::mat zetay = zeta.subvec(m, m + k - 1);
  const arma::mat kronXi = arma::kron(M->Ie, xi);

  const arma::vec eta = arma::solve(
      M->Ie - M->GammaEta - kronXi.t() * M->OmegaEtaXi,
      M->alpha + M->GammaXi * xi + kronXi.t() * M->OmegaXiXi * xi + zetay
  );

  return arma::join_cols(xi, eta);
}


// [[Rcpp::export]]
arma::vec impliedYFromEtaCpp(
    const arma::vec eta,
    const SEXP xptr
) {
  const Rcpp::XPtr<ModelMatrices> M(xptr);
  return M->tau + M->Lambda * eta;
}


// [[Rcpp::export]]
double logLikFromZetaMLCpp(
    const arma::vec zeta,
    const arma::vec y,
    const SEXP xptr,
    const arma::uvec idx
) {
  const Rcpp::XPtr<ModelMatrices> M(xptr);
  const arma::vec eta = impliedEtaFromZetaCpp(zeta, xptr);
  const arma::vec yhat = impliedYFromEtaCpp(eta, xptr);
  const arma::vec epsilon = yhat.elem(idx) - y;
  const arma::vec mu(epsilon.n_elem, arma::fill::zeros);
  const arma::mat Theta = M->Theta.submat(idx, idx);

  return arma::as_scalar(dmvnfast(
    epsilon.t(), mu, Theta, true, 1, false
  ));
}


// [[Rcpp::export]]
double logLikFromZetaEBMCpp(
  const arma::vec zeta,
  const arma::vec y,
  const SEXP xptr,
  const arma::uvec idx
  ) {
  const Rcpp::XPtr<ModelMatrices> M(xptr);
  const arma::vec mu(zeta.n_elem, arma::fill::zeros);

  const double prior = arma::as_scalar(dmvnfast(zeta.t(), mu, M->Psi, true, 1, false));
  const double evidence = logLikFromZetaMLCpp(zeta, y, xptr, idx);

  return prior + evidence;
}
