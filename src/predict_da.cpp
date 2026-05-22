#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat predictedLatentScoresCpp(
    const arma::mat Zeta,
    const Rcpp::List matrices
) {
  
  // model matrices
  const arma::mat GammaXi    = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat GammaEta   = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat OmegaXiXi  = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat OmegaEtaXi = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  const arma::vec beta0      = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha      = Rcpp::as<arma::vec>(matrices["alpha"]);
  const arma::mat Ie         = Rcpp::as<arma::mat>(matrices["Ieta"]);

  // dim
  const int k = GammaEta.n_cols;
  const int m = GammaXi.n_cols;
  const int n = Zeta.n_rows;

  // Split Zeta into exogenous and endogenous blocks
  const arma::mat Zetay = Zeta.cols(m, m + k - 1L);

  arma::mat Xi = Zeta.cols(0, m - 1L);
  arma::mat Eta = arma::zeros<arma::mat>(n, k);

  Xi.each_row() += beta0.t();

  for (int i = 0; i < n; i++) {

    const arma::vec zetay_i = Zetay.row(i).t();
    const arma::vec xi_i = Xi.row(i).t();
    const arma::mat kronXi_i = arma::kron(Ie, xi_i);
    
    Eta.row(i) += arma::solve(
      Ie - GammaEta - kronXi_i.t() * OmegaEtaXi,
      alpha + GammaXi * xi_i + kronXi_i.t() * OmegaXiXi * xi_i + zetay_i
    ).t();

  }

  return arma::join_rows(Xi, Eta);
}
