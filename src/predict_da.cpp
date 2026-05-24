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
arma::vec gradLogLikFromZetaMLCpp(
  const arma::vec zeta,
  const arma::vec y,
  const SEXP xptr,
  const arma::uvec idx,
  const double eps = 1e-4
) {
  // avoid overhead of getting numerical gradient from within R
  const double ll0 = logLikFromZetaMLCpp(zeta, y, xptr, idx);

  arma::vec grad(zeta.n_elem);
  for (int i = 0; i < (int)grad.n_elem; i++) {
    arma::vec zetai = zeta;
    zetai(i) = zetai(i) + eps;
    grad(i) = (logLikFromZetaMLCpp(zetai, y, xptr, idx) - ll0) / eps;
  }

  return grad;
}


template<class F>
Rcpp::List hessLogLikFromZetaQuadraticFit(
    const arma::vec& zeta,
    F&& fun,
    const double relStep,
    const double minAbsPar
) {
  const std::size_t p = zeta.n_elem;

  if (p == 0) {
    return Rcpp::List::create(
      Rcpp::Named("mean")     = fun(zeta),
      Rcpp::Named("gradient") = arma::vec(),
      Rcpp::Named("Hessian")  = arma::mat(0, 0)
    );
  }

  arma::vec incr = arma::max(arma::abs(zeta), arma::vec(p).fill(minAbsPar)) * relStep;
  for (std::size_t i = 0; i < p; ++i)
    if (incr[i] == 0.0) incr[i] = relStep;

  std::vector<arma::vec> disp;
  disp.reserve(1 + 2 * p + (p * (p - 1)) / 2);
  disp.emplace_back(arma::zeros<arma::vec>(p));

  for (std::size_t i = 0; i < p; ++i) {
    arma::vec v = arma::zeros<arma::vec>(p);
    v[i] =  1.0; disp.push_back(v);
    v[i] = -1.0; disp.push_back(v);
  }

  for (std::size_t i = 0; i < p - 1; ++i) {
    for (std::size_t j = i + 1; j < p; ++j) {
      arma::vec v = arma::zeros<arma::vec>(p);
      v[i] = 1.0;
      v[j] = 1.0;
      disp.push_back(v);
    }
  }

  const std::size_t m = disp.size();
  arma::vec values(m);
  for (std::size_t k = 0; k < m; ++k)
    values[k] = fun(zeta + disp[k] % incr);

  const std::size_t q = 1 + 2 * p + (p * (p - 1)) / 2;
  arma::mat X(m, q, arma::fill::ones);

  std::size_t colId = 1;
  for (std::size_t j = 0; j < p; ++j, ++colId)
    for (std::size_t k = 0; k < m; ++k)
      X(k, colId) = disp[k][j];

  for (std::size_t j = 0; j < p; ++j, ++colId)
    for (std::size_t k = 0; k < m; ++k)
      X(k, colId) = disp[k][j] * disp[k][j];

  for (std::size_t i = 0; i < p - 1; ++i)
    for (std::size_t j = i + 1; j < p; ++j, ++colId)
      for (std::size_t k = 0; k < m; ++k)
        X(k, colId) = disp[k][i] * disp[k][j];

  arma::vec frac(q, arma::fill::ones);
  for (std::size_t j = 0; j < p; ++j)
    frac[1 + j] = incr[j];

  for (std::size_t j = 0; j < p; ++j)
    frac[1 + p + j] = incr[j] * incr[j];

  colId = 1 + 2 * p;
  for (std::size_t i = 0; i < p - 1; ++i)
    for (std::size_t j = i + 1; j < p; ++j, ++colId)
      frac[colId] = incr[i] * incr[j];

  arma::vec coef = arma::solve(X, values) / frac;

  arma::vec grad = coef.subvec(1, p);
  arma::mat Hess(p, p, arma::fill::zeros);

  for (std::size_t j = 0; j < p; ++j)
    Hess(j, j) = 2.0 * coef[1 + p + j];

  colId = 1 + 2 * p;
  for (std::size_t i = 0; i < p - 1; ++i) {
    for (std::size_t j = i + 1; j < p; ++j, ++colId) {
      Hess(i, j) = coef[colId];
      Hess(j, i) = coef[colId];
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mean")     = coef[0],
    Rcpp::Named("gradient") = grad,
    Rcpp::Named("Hessian")  = Hess
  );
}


template<class F>
Rcpp::List hessLogLikFromZetaCpp(
    const arma::vec& zeta,
    F&& fun,
    const double relStep,
    const double minAbsPar
) {
  return hessLogLikFromZetaQuadraticFit(zeta, std::forward<F>(fun), relStep, minAbsPar);
}


// [[Rcpp::export]]
Rcpp::List hessLogLikFromZetaMLCpp(
  const arma::vec zeta,
  const arma::vec y,
  const SEXP xptr,
  const arma::uvec idx,
  const double relStep = 1e-4,
  const double minAbsPar = 0.0
) {
  auto fun = [&](const arma::vec& zi) -> double {
    return logLikFromZetaMLCpp(zi, y, xptr, idx);
  };

  return hessLogLikFromZetaCpp(zeta, fun, relStep, minAbsPar);
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


// [[Rcpp::export]]
arma::vec gradLogLikFromZetaEBMCpp(
  const arma::vec zeta,
  const arma::vec y,
  const SEXP xptr,
  const arma::uvec idx,
  const double eps = 1e-4
) {
  // avoid overhead of getting numerical gradient from within R
  const double ll0 = logLikFromZetaEBMCpp(zeta, y, xptr, idx);

  arma::vec grad(zeta.n_elem);
  for (int i = 0; i < (int)grad.n_elem; i++) {
    arma::vec zetai = zeta;
    zetai(i) = zetai(i) + eps;
    grad(i) = (logLikFromZetaEBMCpp(zetai, y, xptr, idx) - ll0) / eps;
  }

  return grad;
}


// [[Rcpp::export]]
Rcpp::List hessLogLikFromZetaEBMCpp(
  const arma::vec zeta,
  const arma::vec y,
  const SEXP xptr,
  const arma::uvec idx,
  const double relStep = 1e-4,
  const double minAbsPar = 0.0
) {
  auto fun = [&](const arma::vec& zi) -> double {
    return logLikFromZetaEBMCpp(zi, y, xptr, idx);
  };

  return hessLogLikFromZetaCpp(zeta, fun, relStep, minAbsPar);
}
