#include <RcppArmadillo.h>

#include "lms_model.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Implied moments for PML.
//
// PML needs exactly what the LMS kernel already computes: conditional on the
// nonlinear latents the model is linear-Gaussian, so `LMSModel::muSigma(z)`
// gives the mean and covariance of the underlying variables at a node, and a
// pair probability is a rectangle under the corresponding bivariate marginal.
//
// This exists so there is ONE definition of those moments. The R
// reimplementation it replaces had already drifted -- it conditioned on xi
// values while being fed standard normal innovations, ignored omegaEtaXi, and
// only ever worked for a single eta. Routing through LMSModel means PML inherits
// every fix the LMS kernel gets, composites included.


// Moments at each row of `nodes`, one list element per node.
// [[Rcpp::export]]
Rcpp::List pmlMomentsCpp(const Rcpp::List& modFilled, const arma::mat& nodes) {
  LMSModel model(modFilled);

  Rcpp::List out(nodes.n_rows);
  for (arma::uword i = 0; i < nodes.n_rows; ++i) {
    const std::pair<arma::vec, arma::mat> ms = model.muSigma(nodes.row(i).t());
    out[i] = Rcpp::List::create(Rcpp::Named("mean") = ms.first,
                                Rcpp::Named("cov")  = ms.second);
  }

  return out;
}


// Moments with NOTHING conditioned on: k = 0 makes Oi the identity and zeroes
// the node vector, which is the marginal over all latents of the model with its
// product terms dropped.
//
// Those product terms do NOT vanish from this call -- kron(Ie, beta0) is
// generally nonzero -- so the affected etas' entries are wrong here. That is
// harmless because only interaction-free pairs read this: affectedness is closed
// under the structural paths, so (I - Ge - kronZ'Oex) is block triangular and
// the clean rows of its inverse carry no weight on any affected eta. The R side
// is what guarantees the restriction; see pmlCleanIndicators().
// [[Rcpp::export]]
Rcpp::List pmlUnconditionalMomentsCpp(const Rcpp::List& modFilled) {
  LMSModel model(modFilled);
  model.k = 0;
  model.updateCache();

  const std::pair<arma::vec, arma::mat> ms = model.muSigma(arma::vec());

  return Rcpp::List::create(Rcpp::Named("mean") = ms.first,
                            Rcpp::Named("cov")  = ms.second);
}
