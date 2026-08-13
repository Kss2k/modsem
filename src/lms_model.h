#ifndef LMS_MODEL_H
#define LMS_MODEL_H

// The LMS model: parameter matrices, the cached z-independent derived
// quantities, and `muSigma(z)` -- the implied mean and covariance of the
// observed variables at a quadrature node.
//
// Extracted verbatim from equations_lms.cpp so estimators other than the
// full-information likelihood can reuse it. PML (src/pml.cpp) needs exactly
// `muSigma(z)`: conditional on the nonlinear latents the model is
// linear-Gaussian, so those moments are all a pair probability requires.
//
// Include AFTER RcppArmadillo.


inline arma::mat makeOi(unsigned k, unsigned numXis) {
  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));

  return Oi;
}


inline arma::vec makeZvec(unsigned k, unsigned numXis, const arma::vec& z) {
  if (k > 0) return arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else       return arma::zeros<arma::vec>(numXis);
}


struct LMSModel {
  arma::mat A, Oxx, Oex, Ie, lY, lX, W, T, tY, tX, Gx, Ge,
    a, beta0, Psi, d, e, covZetaXi;
  unsigned  k        = 0;
  unsigned  numXis   = 0;
  bool hasComposites = false;

  // Precomputed z-independent derived quantities; call updateCache() after any
  // change to A, Gx, Psi, or covZetaXi (i.e. after lmsParam / setParams perturbations).
  arma::mat Oi, GxA, AOi, AOiAt, ZetaProj, PsiOrth;

  void updateCache() {
    Oi       = makeOi(k, numXis);
    ZetaProj = arma::solve(arma::trimatl(A), covZetaXi.t(), arma::solve_opts::fast).t();
    PsiOrth  = Psi - ZetaProj * ZetaProj.t();
    GxA      = Gx * A;
    AOi      = A * Oi;
    AOiAt    = AOi * A.t();
  }

  explicit LMSModel(const Rcpp::List& modFilled) {

    Rcpp::List matrices = modFilled["matrices"];
    Rcpp::List info     = modFilled["info"];
    Rcpp::List quad     = modFilled["quad"];

    k       = Rcpp::as<unsigned>(quad["k"]);
    numXis  = Rcpp::as<unsigned>(info["numXis"]);
    hasComposites = Rcpp::as<bool>(info["hasComposites"]);

    // one-liners, no loops
    A         = Rcpp::as<arma::mat>(matrices["A"]);
    Oxx       = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
    Oex       = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
    Ie        = Rcpp::as<arma::mat>(matrices["Ieta"]);
    lY        = Rcpp::as<arma::mat>(matrices["lambdaY"]);
    lX        = Rcpp::as<arma::mat>(matrices["lambdaX"]);
    tY        = Rcpp::as<arma::mat>(matrices["tauY"]);
    tX        = Rcpp::as<arma::mat>(matrices["tauX"]);
    W         = Rcpp::as<arma::mat>(matrices["W"]);
    T         = Rcpp::as<arma::mat>(matrices["T"]);
    Gx        = Rcpp::as<arma::mat>(matrices["gammaXi"]);
    Ge        = Rcpp::as<arma::mat>(matrices["gammaEta"]);
    a         = Rcpp::as<arma::mat>(matrices["alpha"]);
    beta0     = Rcpp::as<arma::mat>(matrices["beta0"]);
    Psi       = Rcpp::as<arma::mat>(matrices["psi"]);
    d         = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
    e         = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);
    covZetaXi = Rcpp::as<arma::mat>(matrices["covZetaXi"]);

    updateCache();
  }

  // Combined mu+Sigma: avoids recomputing zVec/kronZ/Binv/lXc for the same z.
  // Use this in any hot path that needs both.
  std::pair<arma::vec, arma::mat> muSigma(const arma::vec& z) const {
    const arma::vec zVec  = makeZvec(k, numXis, z);
    const arma::vec muXi  = beta0 + A * zVec;
    const arma::mat kronZ = arma::kron(Ie, muXi);
    const arma::mat B     = Ie - Ge - kronZ.t() * Oex;

    const arma::mat Binv     = arma::inv(B);
    const arma::vec muEta    = Binv * (a + Gx * muXi + kronZ.t() * Oxx * muXi + ZetaProj * zVec);
    const arma::mat Eta      = Binv * (GxA + kronZ.t() * Oxx * A + ZetaProj);
    const arma::mat varXi    = AOiAt;
    const arma::mat varEta   = Eta * Oi * Eta.t() + Binv * PsiOrth * Binv.t();
    const arma::mat covXiEta = AOi * Eta.t();

    const arma::mat vcovXiEta = arma::join_cols(
      arma::join_rows(varXi,        covXiEta),
      arma::join_rows(covXiEta.t(), varEta)
    );

    const arma::vec xieta = arma::join_cols(muXi, muEta);

    if (hasComposites) {
      const arma::mat lXc = lX + T * W * arma::pinv(W.t() * T * W);
      const arma::mat dc  = d + T - lXc * W.t() * T * W * lXc.t();
      return std::make_pair(tX + lXc * xieta, lXc * vcovXiEta * lXc.t() + dc);
    }

    return std::make_pair(tX + lX * xieta, lX * vcovXiEta * lX.t() + d);
  }

  LMSModel threadClone() const {
    LMSModel c = *this;    // shallow for everything (fast)
                           // Deep-copy ONLY what setParams()/lmsParam can modify:
    c.A         = arma::mat(A);
    c.Oxx       = arma::mat(Oxx);
    c.Oex       = arma::mat(Oex);
    c.Ie        = arma::mat(Ie);
    c.lY        = arma::mat(lY);
    c.lX        = arma::mat(lX);
    c.tY        = arma::mat(tY);
    c.tX        = arma::mat(tX);
    c.W         = arma::mat(W);
    c.T         = arma::mat(T);
    c.Gx        = arma::mat(Gx);
    c.Ge        = arma::mat(Ge);
    c.a         = arma::mat(a);
    c.beta0     = arma::mat(beta0);
    c.Psi       = arma::mat(Psi);
    c.covZetaXi = arma::mat(covZetaXi);
    c.d         = arma::mat(d);
    c.e         = arma::mat(e);
    c.Oi        = arma::mat(Oi);
    c.GxA       = arma::mat(GxA);
    c.AOi       = arma::mat(AOi);
    c.AOiAt     = arma::mat(AOiAt);
    c.ZetaProj  = arma::mat(ZetaProj);
    c.PsiOrth   = arma::mat(PsiOrth);

    return c;
  }
};
#endif
