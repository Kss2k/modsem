#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>

#include "lms.h"
#include "utils.h"
#include "mvnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]


// Deprecated will remove soon...
// [[Rcpp::export]]
arma::vec muLmsCpp(Rcpp::List model, arma::vec z) {
  const Rcpp::List matrices = model["matrices"];
  const Rcpp::List info = model["info"];
  const Rcpp::List quad = model["quad"];
  const int numXis = Rcpp::as<int>(info["numXis"]);
  const int k = Rcpp::as<int>(quad["k"]);
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  const arma::mat Ie = Rcpp::as<arma::mat>(matrices["Ieta"]);
  const arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  const arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  const arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else       zVec = arma::zeros<arma::vec>(numXis);

  const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);
  const arma::mat Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

  const arma::vec muXi  = beta0 + A * zVec;
  const arma::vec muEta = Binv * (a +
      Gx * (beta0 + A * zVec) +
      kronZ.t() * Oxx * (beta0 + A * zVec));

  return tX + lX * arma::join_cols(muXi, muEta);
}


// Deprecated will remove soon...
// [[Rcpp::export]]
arma::mat sigmaLmsCpp(Rcpp::List model, arma::vec z) {
  const Rcpp::List matrices = model["matrices"];
  const Rcpp::List info = model["info"];
  const Rcpp::List quad = model["quad"];
  const int numXis = Rcpp::as<int>(info["numXis"]);
  const int k = Rcpp::as<int>(quad["k"]);
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  const arma::mat Ie = Rcpp::as<arma::mat>(matrices["Ieta"]);
  const arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  const arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  const arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);
  const arma::mat Psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat d = Rcpp::as<arma::mat>(matrices["thetaDelta"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else       zVec = arma::zeros<arma::vec>(numXis);

  const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);
  const arma::mat Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));

  const arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
  const arma::mat varXi    = A * Oi * A.t();
  const arma::mat varEta   = Eta * Oi * Eta.t() + Binv * Psi * Binv.t();
  const arma::mat covXiEta = A * Oi * Eta.t();

  const arma::mat vcovXiEta = arma::join_cols(
    arma::join_rows(varXi, covXiEta),
    arma::join_rows(covXiEta.t(), varEta)
  );

  return lX * vcovXiEta * lX.t() + d;
}


inline arma::mat make_Oi(unsigned k, unsigned numXis) {
  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));

  return Oi;
}


inline arma::vec make_zvec(unsigned k, unsigned numXis, const arma::vec& z) {
  if (k > 0) return arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else       return arma::zeros<arma::vec>(numXis);
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
    const arma::vec zVec = make_zvec(k, numXis, z);
    const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);
    const arma::mat Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    const arma::vec muXi  = beta0 + A * zVec;
    const arma::vec muEta = Binv * (a +
        Gx * (beta0 + A * zVec) +
        kronZ.t() * Oxx * (beta0 + A * zVec));

    return tX + lX * arma::join_cols(muXi, muEta);
  }

  arma::mat Sigma(const arma::vec& z) const {
    const arma::vec zVec  = make_zvec(k, numXis, z);
    const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);
    const arma::mat Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    const arma::mat Oi = make_Oi(k, numXis);
    const arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
    const arma::mat varXi    = A * Oi * A.t();
    const arma::mat varEta   = Eta * Oi * Eta.t() + Binv * Psi * Binv.t();
    const arma::mat covXiEta = A * Oi * Eta.t();

    const arma::mat vcovXiEta = arma::join_cols(
        arma::join_rows(varXi, covXiEta),
        arma::join_rows(covXiEta.t(), varEta)
        );

    return lX * vcovXiEta * lX.t() + d;
  }

  LMSModel thread_clone() const {
    LMSModel c = *this;    // shallow for everything (fast)
                           // Deep-copy ONLY what setParams()/lms_param can modify:
    c.A     = arma::mat(A);
    c.Oxx   = arma::mat(Oxx);
    c.Oex   = arma::mat(Oex);
    c.Ie    = arma::mat(Ie);
    c.lY    = arma::mat(lY);
    c.lX    = arma::mat(lX);
    c.tY    = arma::mat(tY);
    c.tX    = arma::mat(tX);
    c.Gx    = arma::mat(Gx);
    c.Ge    = arma::mat(Ge);
    c.a     = arma::mat(a);
    c.beta0 = arma::mat(beta0);
    c.Psi   = arma::mat(Psi);
    c.d     = arma::mat(d);
    c.e     = arma::mat(e);

    return c;
  }
};

inline double& lms_param(LMSModel& M, std::size_t blk,
          std::size_t r, std::size_t c) {
  switch (blk) {
    case 0 : return  M.lX   (r,c);
    case 1 : return  M.lY   (r,c);
    case 2 : return  M.tX   (r,c);
    case 3 : return  M.tY   (r,c);
    case 4 : return  M.d    (r,c);
    case 5 : return  M.e    (r,c);
    case 6 : return  M.A    (r,c);
    case 7 : return  M.Psi  (r,c);
    case 8 : return  M.a    (r,c);
    case 9 : return  M.beta0(r,c);
    case 10: return  M.Gx   (r,c);
    case 11: return  M.Ge   (r,c);
    case 12: return  M.Oxx  (r,c);
    case 13: return  M.Oex  (r,c);
    default: Rcpp::stop("unknown block id");
  }
}


arma::mat gradAlphaComplete(const LMSModel&  M,
                            const arma::mat& V,
                            const std::vector<arma::vec>& TGamma,
                            const std::vector<std::vector<arma::vec>>& MeanPatterns,
                            const std::vector<std::vector<arma::mat>>& CovPatterns,
                            const std::vector<arma::uvec>& colidx,
                            const int npatterns);

arma::mat gradBeta0Complete(const LMSModel&  M,
                            const arma::mat& V,
                            const std::vector<arma::vec>& TGamma,
                            const std::vector<std::vector<arma::vec>>& MeanPatterns,
                            const std::vector<std::vector<arma::mat>>& CovPatterns,
                            const std::vector<arma::uvec>& colidx,
                            const int npatterns);

arma::mat gradTauXComplete(const LMSModel&  M,
                           const arma::mat& V,
                           const std::vector<arma::vec>& TGamma,
                           const std::vector<std::vector<arma::vec>>& MeanPatterns,
                           const std::vector<std::vector<arma::mat>>& CovPatterns,
                           const std::vector<arma::uvec>& colidx,
                           const int npatterns);

arma::mat gradLambdaXComplete(const LMSModel&  M,
                              const arma::mat& V,
                              const std::vector<arma::vec>& TGamma,
                              const std::vector<std::vector<arma::vec>>& MeanPatterns,
                              const std::vector<std::vector<arma::mat>>& CovPatterns,
                              const std::vector<arma::uvec>& colidx,
                              const int npatterns);

arma::mat gradAComplete(const LMSModel&  M,
                        const arma::mat& V,
                        const std::vector<arma::vec>& TGamma,
                        const std::vector<std::vector<arma::vec>>& MeanPatterns,
                        const std::vector<std::vector<arma::mat>>& CovPatterns,
                        const std::vector<arma::uvec>& colidx,
                        const int npatterns);

arma::mat gradThetaDeltaComplete(const LMSModel&  M,
                                 const arma::mat& V,
                                 const std::vector<arma::vec>& TGamma,
                                 const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                 const std::vector<std::vector<arma::mat>>& CovPatterns,
                                 const std::vector<arma::uvec>& colidx,
                                 const int npatterns);

arma::mat gradGammaXiComplete(const LMSModel&  M,
                              const arma::mat& V,
                              const std::vector<arma::vec>& TGamma,
                              const std::vector<std::vector<arma::vec>>& MeanPatterns,
                              const std::vector<std::vector<arma::mat>>& CovPatterns,
                              const std::vector<arma::uvec>& colidx,
                              const int npatterns);

arma::mat gradGammaEtaComplete(const LMSModel&  M,
                               const arma::mat& V,
                               const std::vector<arma::vec>& TGamma,
                               const std::vector<std::vector<arma::vec>>& MeanPatterns,
                               const std::vector<std::vector<arma::mat>>& CovPatterns,
                               const std::vector<arma::uvec>& colidx,
                               const int npatterns);

arma::mat gradOmegaXiXiComplete(const LMSModel&  M,
                                const arma::mat& V,
                                const std::vector<arma::vec>& TGamma,
                                const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                const std::vector<std::vector<arma::mat>>& CovPatterns,
                                const std::vector<arma::uvec>& colidx,
                                const int npatterns);

arma::mat gradOmegaEtaXiComplete(const LMSModel&  M,
                                 const arma::mat& V,
                                 const std::vector<arma::vec>& TGamma,
                                 const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                 const std::vector<std::vector<arma::mat>>& CovPatterns,
                                 const std::vector<arma::uvec>& colidx,
                                 const int npatterns);

arma::mat gradLambdaXComplete(const LMSModel&  M,
                              const arma::mat& V,
                              const std::vector<arma::vec>& TGamma,
                              const std::vector<std::vector<arma::vec>>& MeanPatterns,
                              const std::vector<std::vector<arma::mat>>& CovPatterns,
                              const std::vector<arma::uvec>& colidx,
                              const int npatterns) {
  arma::mat gradLambda(M.lX.n_rows, M.lX.n_cols, arma::fill::zeros);

  const std::size_t J = V.n_rows;
  const arma::mat Oi = make_Oi(M.k, M.numXis);

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv = arma::inv(B);
    const arma::vec rhs = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;
    const arma::vec latent = arma::join_cols(muXi, muEta);
    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);
    const arma::mat Eta = Binv * (M.Gx * M.A + kronZ.t() * M.Oxx * M.A);
    const arma::mat varXi = M.A * Oi * M.A.t();
    const arma::mat varEta = Eta * Oi * Eta.t() + Binv * M.Psi * Binv.t();
    const arma::mat covXiEta = M.A * Oi * Eta.t();

    arma::mat vcov(varXi.n_rows + varEta.n_rows,
                   varXi.n_cols + varEta.n_cols, arma::fill::zeros);
    vcov.submat(0, 0, varXi.n_rows - 1, varXi.n_cols - 1) = varXi;
    vcov.submat(0, varXi.n_cols, varXi.n_rows - 1, vcov.n_cols - 1) = covXiEta;
    vcov.submat(varXi.n_rows, 0, vcov.n_rows - 1, varXi.n_cols - 1) = covXiEta.t();
    vcov.submat(varXi.n_rows, varXi.n_cols, vcov.n_rows - 1, vcov.n_cols - 1) = varEta;

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);
    std::vector<arma::mat> gradSigmaPatterns(npatterns);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) {
        gradSigmaPatterns[i].reset();
        continue;
      }

      const arma::uvec& idx = colidx[i];
      const arma::vec& nu = MeanPatterns[j][i];
      const arma::mat& S  = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) {
        gradSigmaPatterns[i].reset();
        continue;
      }

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);

      gradMu.elem(idx) += tg * linv_diff;

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();
      arma::mat termS = SigInv * S * SigInv;
      arma::mat termDiff = SigInv * diffOuter * SigInv;
      gradSigmaPatterns[i] = 0.5 * (termS + tg * termDiff - tg * SigInv);
    }

    if (arma::all(gradMu == 0.0) &&
        std::all_of(gradSigmaPatterns.begin(), gradSigmaPatterns.end(),
                    [](const arma::mat& m){ return m.is_empty(); })) continue;

    for (arma::uword r = 0; r < M.lX.n_rows; ++r) {
      for (arma::uword c = 0; c < M.lX.n_cols; ++c) {
        double mean_contrib = gradMu[r] * latent[c];

        arma::mat E(M.lX.n_rows, M.lX.n_cols, arma::fill::zeros);
        E(r, c) = 1.0;
        arma::mat dSigma = E * vcov * M.lX.t() + M.lX * vcov * E.t();

        double cov_contrib = 0.0;
        for (int i = 0; i < npatterns; ++i) {
          if (gradSigmaPatterns[i].is_empty()) continue;
          const arma::uvec& idx = colidx[i];
          if (idx.is_empty()) continue;
          arma::mat dSigmaSub = dSigma.submat(idx, idx);
          cov_contrib += arma::accu(gradSigmaPatterns[i] % dSigmaSub);
        }

        gradLambda(r, c) += mean_contrib + cov_contrib;
      }
    }
  }

  return gradLambda;
}

arma::mat gradAComplete(const LMSModel&  M,
                        const arma::mat& V,
                        const std::vector<arma::vec>& TGamma,
                        const std::vector<std::vector<arma::vec>>& MeanPatterns,
                        const std::vector<std::vector<arma::mat>>& CovPatterns,
                        const std::vector<arma::uvec>& colidx,
                        const int npatterns) {
  arma::mat gradA(M.A.n_rows, M.A.n_cols, arma::fill::zeros);
  const arma::mat Oi = make_Oi(M.k, M.numXis);

  const std::size_t J = V.n_rows;
  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B     = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv  = arma::inv(B);
    const arma::vec rhs   = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;

    const arma::vec latent = arma::join_cols(muXi, muEta);
    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);
    arma::mat gradSigma(Sigma.n_rows, Sigma.n_cols, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;
      const arma::uvec& idx = colidx[i];
      const arma::vec& nu  = MeanPatterns[j][i];
      const arma::mat& S   = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);
      gradMu.elem(idx) += tg * linv_diff;

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();
      arma::mat gradSel = 0.5 * (SigInv * S * SigInv
                                 + tg * SigInv * diffOuter * SigInv
                                 - tg * SigInv);
      gradSigma.submat(idx, idx) += gradSel;
    }

    const arma::mat Q   = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat Eta = Binv * Q;
    const arma::uword numXi  = muXi.n_elem;
    const arma::uword numEta = muEta.n_elem;

    for (arma::uword r = 0; r < M.A.n_rows; ++r) {
      arma::vec e_r(M.A.n_rows, arma::fill::zeros); e_r[r] = 1.0;
      for (arma::uword c = 0; c < M.A.n_cols; ++c) {
        arma::vec e_c(M.A.n_cols, arma::fill::zeros); e_c[c] = 1.0;
        arma::mat dA = e_r * e_c.t();

        arma::vec d_muXi(numXi, arma::fill::zeros);
        d_muXi[r] = zVec[c];
        arma::mat d_kronZ = arma::kron(M.Ie, d_muXi);
        arma::mat dB = -d_kronZ.t() * M.Oex;

        arma::vec d_rhs = M.Gx * d_muXi +
                          d_kronZ.t() * M.Oxx * muXi +
                          kronZ.t() * M.Oxx * d_muXi;
        arma::vec d_muEta = Binv * (d_rhs - dB * muEta);
        arma::mat dBinv = -Binv * dB * Binv;

        arma::mat dVarXi = dA * Oi * M.A.t() + M.A * Oi * dA.t();
        arma::mat dQ = M.Gx * dA + d_kronZ.t() * M.Oxx * M.A + kronZ.t() * M.Oxx * dA;
        arma::mat dEta = dBinv * Q + Binv * dQ;

        arma::mat dVarEta = dEta * Oi * Eta.t() + Eta * Oi * dEta.t() +
                            dBinv * M.Psi * Binv.t() + Binv * M.Psi * dBinv.t();
        arma::mat dCovXiEta = dA * Oi * Eta.t() + M.A * Oi * dEta.t();

        arma::mat dVcov(numXi + numEta, numXi + numEta, arma::fill::zeros);
        dVcov.submat(0, 0, numXi - 1, numXi - 1) = dVarXi;
        dVcov.submat(0, numXi, numXi - 1, numXi + numEta - 1) = dCovXiEta;
        dVcov.submat(numXi, 0, numXi + numEta - 1, numXi - 1) = dCovXiEta.t();
        dVcov.submat(numXi, numXi, numXi + numEta - 1, numXi + numEta - 1) = dVarEta;

        arma::vec dLatent = arma::join_cols(d_muXi, d_muEta);
        arma::vec dMu = M.lX * dLatent;
        arma::mat dSigma = M.lX * dVcov * M.lX.t();

        const double mean_contrib = arma::dot(gradMu, dMu);
        const double cov_contrib  = arma::accu(gradSigma % dSigma);
        gradA(r, c) += mean_contrib + cov_contrib;
      }
    }
  }

  return gradA;
}

arma::mat gradPsiComplete(const LMSModel&  M,
                          const arma::mat& V,
                          const std::vector<arma::vec>& TGamma,
                          const std::vector<std::vector<arma::vec>>& MeanPatterns,
                          const std::vector<std::vector<arma::mat>>& CovPatterns,
                          const std::vector<arma::uvec>& colidx,
                          const int npatterns) {
  arma::mat gradPsi(M.Psi.n_rows, M.Psi.n_cols, arma::fill::zeros);

  const std::size_t J = V.n_rows;
  const arma::uword dimXi  = M.numXis;
  const arma::uword dimEta = M.Ie.n_rows;

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv = arma::inv(B);

    const arma::vec mu = M.mu(z);
    const arma::mat Sigma = M.Sigma(z);

    arma::mat gradSigma(Sigma.n_rows, Sigma.n_cols, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;

      const arma::uvec& idx = colidx[i];
      const arma::vec& nu  = MeanPatterns[j][i];
      const arma::mat& S   = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();

      arma::mat gradSel = 0.5 * (SigInv * S * SigInv
                                 + tg * SigInv * diffOuter * SigInv
                                 - tg * SigInv);

      gradSigma.submat(idx, idx) += gradSel;
    }

    if (!gradSigma.is_finite()) continue;
    if (arma::approx_equal(gradSigma, arma::zeros<arma::mat>(gradSigma.n_rows, gradSigma.n_cols), "absdiff", 0.0))
      continue;

    const arma::mat Gvcov = M.lX.t() * gradSigma * M.lX;
    const arma::mat GvarEta = Gvcov.submat(dimXi, dimXi,
                                          dimXi + dimEta - 1,
                                          dimXi + dimEta - 1);

    gradPsi += Binv.t() * GvarEta * Binv;
  }

  return gradPsi;
}


arma::mat gradThetaDeltaComplete(const LMSModel&  M,
                                 const arma::mat& V,
                                 const std::vector<arma::vec>& TGamma,
                                 const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                 const std::vector<std::vector<arma::mat>>& CovPatterns,
                                 const std::vector<arma::uvec>& colidx,
                                 const int npatterns) {
  arma::mat gradTheta = arma::zeros<arma::mat>(M.d.n_rows, M.d.n_cols);

  const std::size_t J = V.n_rows;

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z  = V.row(j).t();
    const arma::vec mu = M.mu(z);
    const arma::mat Sigma = M.Sigma(z);

    arma::mat gradSigma(Sigma.n_rows, Sigma.n_cols, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;

      const arma::uvec& idx = colidx[i];
      const arma::vec& nu  = MeanPatterns[j][i];
      const arma::mat& S   = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();

      arma::mat gradSel = 0.5 * (SigInv * S * SigInv
                                 + tg * SigInv * diffOuter * SigInv
                                 - tg * SigInv);

      gradSigma.submat(idx, idx) += gradSel;
    }

    gradTheta += gradSigma;
  }

  return gradTheta;
}


arma::mat gradGammaXiComplete(const LMSModel&  M,
                              const arma::mat& V,
                              const std::vector<arma::vec>& TGamma,
                              const std::vector<std::vector<arma::vec>>& MeanPatterns,
                              const std::vector<std::vector<arma::mat>>& CovPatterns,
                              const std::vector<arma::uvec>& colidx,
                              const int npatterns) {
  arma::mat gradGx(M.Gx.n_rows, M.Gx.n_cols, arma::fill::zeros);
  const arma::mat Oi = make_Oi(M.k, M.numXis);

  const std::size_t J = V.n_rows;
  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B     = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv  = arma::inv(B);
    const arma::vec rhs   = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;
    const arma::vec latent = arma::join_cols(muXi, muEta);

    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);
    arma::mat gradSigma(Sigma.n_rows, Sigma.n_cols, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;
      const arma::uvec& idx = colidx[i];
      const arma::vec& nu  = MeanPatterns[j][i];
      const arma::mat& S   = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);
      gradMu.elem(idx) += tg * linv_diff;

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();
      arma::mat gradSel = 0.5 * (SigInv * S * SigInv
                                 + tg * SigInv * diffOuter * SigInv
                                 - tg * SigInv);
      gradSigma.submat(idx, idx) += gradSel;
    }

    const arma::mat Q   = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat Eta = Binv * Q;
    const arma::uword numXi  = muXi.n_elem;
    const arma::uword numEta = muEta.n_elem;

    for (arma::uword r = 0; r < M.Gx.n_rows; ++r) {
      arma::vec e_r(M.Gx.n_rows, arma::fill::zeros); e_r[r] = 1.0;
      for (arma::uword c = 0; c < M.Gx.n_cols; ++c) {
        arma::vec e_c(M.Gx.n_cols, arma::fill::zeros); e_c[c] = 1.0;
        arma::mat dG = e_r * e_c.t();

        arma::vec d_muXi(numXi, arma::fill::zeros);
        arma::vec d_muEta = Binv * (dG * muXi);

        arma::vec dLatent = arma::join_cols(d_muXi, d_muEta);
        arma::vec dMu = M.lX * dLatent;

        arma::mat dVarXi(numXi, numXi, arma::fill::zeros);
        arma::mat dQ = dG * M.A;
        arma::mat dEta = Binv * dQ;
        arma::mat dVarEta = dEta * Oi * Eta.t() + Eta * Oi * dEta.t();
        arma::mat dCovXiEta = M.A * Oi * dEta.t();

        arma::mat dVcov(numXi + numEta, numXi + numEta, arma::fill::zeros);
        dVcov.submat(0, 0, numXi - 1, numXi - 1) = dVarXi;
        dVcov.submat(0, numXi, numXi - 1, numXi + numEta - 1) = dCovXiEta;
        dVcov.submat(numXi, 0, numXi + numEta - 1, numXi - 1) = dCovXiEta.t();
        dVcov.submat(numXi, numXi, numXi + numEta - 1, numXi + numEta - 1) = dVarEta;

        arma::mat dSigma = M.lX * dVcov * M.lX.t();
        const double mean_contrib = arma::dot(gradMu, dMu);
        const double cov_contrib  = arma::accu(gradSigma % dSigma);
        gradGx(r, c) += mean_contrib + cov_contrib;
      }
    }
  }

  return gradGx;
}


arma::mat gradGammaEtaComplete(const LMSModel&  M,
                               const arma::mat& V,
                               const std::vector<arma::vec>& TGamma,
                               const std::vector<std::vector<arma::vec>>& MeanPatterns,
                               const std::vector<std::vector<arma::mat>>& CovPatterns,
                               const std::vector<arma::uvec>& colidx,
                               const int npatterns) {
  arma::mat gradGe(M.Ge.n_rows, M.Ge.n_cols, arma::fill::zeros);
  const arma::mat Oi = make_Oi(M.k, M.numXis);

  const std::size_t J = V.n_rows;
  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B     = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv  = arma::inv(B);
    const arma::vec rhs   = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;
    const arma::vec latent = arma::join_cols(muXi, muEta);

    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);
    arma::mat gradSigma(Sigma.n_rows, Sigma.n_cols, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;
      const arma::uvec& idx = colidx[i];
      const arma::vec& nu  = MeanPatterns[j][i];
      const arma::mat& S   = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);
      gradMu.elem(idx) += tg * linv_diff;

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();
      arma::mat gradSel = 0.5 * (SigInv * S * SigInv
                                 + tg * SigInv * diffOuter * SigInv
                                 - tg * SigInv);
      gradSigma.submat(idx, idx) += gradSel;
    }

    const arma::mat Q   = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat Eta = Binv * Q;
    const arma::uword numXi  = muXi.n_elem;
    const arma::uword numEta = muEta.n_elem;

    for (arma::uword r = 0; r < M.Ge.n_rows; ++r) {
      arma::vec e_r(M.Ge.n_rows, arma::fill::zeros); e_r[r] = 1.0;
      for (arma::uword c = 0; c < M.Ge.n_cols; ++c) {
        arma::vec e_c(M.Ge.n_cols, arma::fill::zeros); e_c[c] = 1.0;
        arma::mat dGe = e_r * e_c.t();

        arma::mat dB = -dGe;
        arma::mat dBinv = Binv * dGe * Binv;
        arma::vec d_muXi(numXi, arma::fill::zeros);
        arma::vec d_muEta = dBinv * rhs;

        arma::vec dLatent = arma::join_cols(d_muXi, d_muEta);
        arma::vec dMu = M.lX * dLatent;

        arma::mat dVarXi(numXi, numXi, arma::fill::zeros);
        arma::mat dEta = dBinv * Q;
        arma::mat dVarEta = dEta * Oi * Eta.t() + Eta * Oi * dEta.t() +
                            dBinv * M.Psi * Binv.t() + Binv * M.Psi * dBinv.t();
        arma::mat dCovXiEta = M.A * Oi * dEta.t();

        arma::mat dVcov(numXi + numEta, numXi + numEta, arma::fill::zeros);
        dVcov.submat(0, 0, numXi - 1, numXi - 1) = dVarXi;
        dVcov.submat(0, numXi, numXi - 1, numXi + numEta - 1) = dCovXiEta;
        dVcov.submat(numXi, 0, numXi + numEta - 1, numXi - 1) = dCovXiEta.t();
        dVcov.submat(numXi, numXi, numXi + numEta - 1, numXi + numEta - 1) = dVarEta;

        arma::mat dSigma = M.lX * dVcov * M.lX.t();
        const double mean_contrib = arma::dot(gradMu, dMu);
        const double cov_contrib  = arma::accu(gradSigma % dSigma);
        gradGe(r, c) += mean_contrib + cov_contrib;
      }
    }
  }

  return gradGe;
}


arma::mat gradOmegaXiXiComplete(const LMSModel&  M,
                                const arma::mat& V,
                                const std::vector<arma::vec>& TGamma,
                                const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                const std::vector<std::vector<arma::mat>>& CovPatterns,
                                const std::vector<arma::uvec>& colidx,
                                const int npatterns) {
  arma::mat gradOxx(M.Oxx.n_rows, M.Oxx.n_cols, arma::fill::zeros);
  const arma::mat Oi = make_Oi(M.k, M.numXis);

  const std::size_t J = V.n_rows;
  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B     = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv  = arma::inv(B);
    const arma::vec rhs   = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;
    const arma::vec latent = arma::join_cols(muXi, muEta);

    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);
    arma::mat gradSigma(Sigma.n_rows, Sigma.n_cols, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;
      const arma::uvec& idx = colidx[i];
      const arma::vec& nu  = MeanPatterns[j][i];
      const arma::mat& S   = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);
      gradMu.elem(idx) += tg * linv_diff;

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();
      arma::mat gradSel = 0.5 * (SigInv * S * SigInv
                                 + tg * SigInv * diffOuter * SigInv
                                 - tg * SigInv);
      gradSigma.submat(idx, idx) += gradSel;
    }

    const arma::mat Q   = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat Eta = Binv * Q;
    const arma::uword numXi  = muXi.n_elem;
    const arma::uword numEta = muEta.n_elem;

    for (arma::uword r = 0; r < M.Oxx.n_rows; ++r) {
      for (arma::uword c = 0; c < M.Oxx.n_cols; ++c) {
        arma::mat dOxx(M.Oxx.n_rows, M.Oxx.n_cols, arma::fill::zeros);
        dOxx(r, c) = 1.0;

        arma::vec d_muXi(numXi, arma::fill::zeros);
        arma::vec d_muEta = Binv * (kronZ.t() * dOxx * muXi);

        arma::vec dLatent = arma::join_cols(d_muXi, d_muEta);
        arma::vec dMu = M.lX * dLatent;

        arma::mat dVarXi(numXi, numXi, arma::fill::zeros);
        arma::mat dQ = kronZ.t() * dOxx * M.A;
        arma::mat dEta = Binv * dQ;
        arma::mat dVarEta = dEta * Oi * Eta.t() + Eta * Oi * dEta.t();
        arma::mat dCovXiEta = M.A * Oi * dEta.t();

        arma::mat dVcov(numXi + numEta, numXi + numEta, arma::fill::zeros);
        dVcov.submat(0, 0, numXi - 1, numXi - 1) = dVarXi;
        dVcov.submat(0, numXi, numXi - 1, numXi + numEta - 1) = dCovXiEta;
        dVcov.submat(numXi, 0, numXi + numEta - 1, numXi - 1) = dCovXiEta.t();
        dVcov.submat(numXi, numXi, numXi + numEta - 1, numXi + numEta - 1) = dVarEta;

        arma::mat dSigma = M.lX * dVcov * M.lX.t();
        const double mean_contrib = arma::dot(gradMu, dMu);
        const double cov_contrib  = arma::accu(gradSigma % dSigma);
        gradOxx(r, c) += mean_contrib + cov_contrib;
      }
    }
  }

  return gradOxx;
}


arma::mat gradOmegaEtaXiComplete(const LMSModel&  M,
                                 const arma::mat& V,
                                 const std::vector<arma::vec>& TGamma,
                                 const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                 const std::vector<std::vector<arma::mat>>& CovPatterns,
                                 const std::vector<arma::uvec>& colidx,
                                 const int npatterns) {
  arma::mat gradOex(M.Oex.n_rows, M.Oex.n_cols, arma::fill::zeros);
  const arma::mat Oi = make_Oi(M.k, M.numXis);

  const std::size_t J = V.n_rows;
  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B     = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv  = arma::inv(B);
    const arma::vec rhs   = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;
    const arma::vec latent = arma::join_cols(muXi, muEta);

    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);
    arma::mat gradSigma(Sigma.n_rows, Sigma.n_cols, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;
      const arma::uvec& idx = colidx[i];
      const arma::vec& nu  = MeanPatterns[j][i];
      const arma::mat& S   = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);
      gradMu.elem(idx) += tg * linv_diff;

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();
      arma::mat gradSel = 0.5 * (SigInv * S * SigInv
                                 + tg * SigInv * diffOuter * SigInv
                                 - tg * SigInv);
      gradSigma.submat(idx, idx) += gradSel;
    }

    const arma::mat Q   = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat Eta = Binv * Q;
    const arma::uword numXi  = muXi.n_elem;
    const arma::uword numEta = muEta.n_elem;

    for (arma::uword r = 0; r < M.Oex.n_rows; ++r) {
      for (arma::uword c = 0; c < M.Oex.n_cols; ++c) {
        arma::mat dOex(M.Oex.n_rows, M.Oex.n_cols, arma::fill::zeros);
        dOex(r, c) = 1.0;

        arma::mat dB = -kronZ.t() * dOex;
        arma::mat dBinv = -Binv * dB * Binv;
        arma::vec d_muXi(numXi, arma::fill::zeros);
        arma::vec d_muEta = dBinv * rhs;

        arma::vec dLatent = arma::join_cols(d_muXi, d_muEta);
        arma::vec dMu = M.lX * dLatent;

        arma::mat dVarXi(numXi, numXi, arma::fill::zeros);
        arma::mat dEta = dBinv * Q;
        arma::mat dVarEta = dEta * Oi * Eta.t() + Eta * Oi * dEta.t() +
                            dBinv * M.Psi * Binv.t() + Binv * M.Psi * dBinv.t();
        arma::mat dCovXiEta = M.A * Oi * dEta.t();

        arma::mat dVcov(numXi + numEta, numXi + numEta, arma::fill::zeros);
        dVcov.submat(0, 0, numXi - 1, numXi - 1) = dVarXi;
        dVcov.submat(0, numXi, numXi - 1, numXi + numEta - 1) = dCovXiEta;
        dVcov.submat(numXi, 0, numXi + numEta - 1, numXi - 1) = dCovXiEta.t();
        dVcov.submat(numXi, numXi, numXi + numEta - 1, numXi + numEta - 1) = dVarEta;

        arma::mat dSigma = M.lX * dVcov * M.lX.t();
        const double mean_contrib = arma::dot(gradMu, dMu);
        const double cov_contrib  = arma::accu(gradSigma % dSigma);
        gradOex(r, c) += mean_contrib + cov_contrib;
      }
    }
  }

  return gradOex;
}


template<class F>
arma::vec gradientFD(LMSModel&         M,
                     F&&               logLik,
                     const arma::uvec& block,
                     const arma::uvec& row,
                     const arma::uvec& col,
                     const arma::uvec& symmetric,
                     const double      eps = 1e-6,
                     const int         ncores = 1L) {
  ThreadSetter ts(ncores);

  const std::size_t p = block.n_elem;
  arma::vec grad(p);

  // Baseline likelihood on the original (unmodified) model:
  const double f0 = logLik(M);

  // Parallelize over coordinates. Each iteration creates its own model copy.
  // NOTE: We mark logLik firstprivate so each thread gets its own copy of the functor/lambda.
  // We only read from M to construct the thread-local copy, so sharing M is OK.
  #pragma omp parallel for default(none) \
      shared(M, block, row, col, symmetric, eps, grad, f0, p) \
      firstprivate(logLik) \
      schedule(static)
  for (std::size_t k = 0; k < p; ++k) {
    // Thread-local model instance
    LMSModel Mc = M.thread_clone();

    // Access the parameter(s) to perturb in the *local* model:
    double& ti   = lms_param(Mc, block[k], row[k], col[k]);
    double* tj   = nullptr;

    if (symmetric[k] && row[k] != col[k]) {
      tj = &lms_param(Mc, block[k], col[k], row[k]); // symmetric partner
    }

    // Forward finite difference step
    ti += eps;
    if (tj) *tj += eps;

    // Evaluate on the perturbed *local* model
    const double f1 = logLik(Mc);

    // Gradient component
    grad[k] = (f1 - f0) / eps;

    // No need to restore: Mc is thread-local and will be destroyed here.
  }

  return grad;
}


template<class F>
arma::vec gradientCentralFD(LMSModel&         M,
                            F&&               logLik,
                            const arma::uvec& block,
                            const arma::uvec& row,
                            const arma::uvec& col,
                            const arma::uvec& symmetric,
                            const double      eps = 1e-6,
                            const int         ncores = 1L) {
  ThreadSetter ts(ncores);

  const std::size_t p = block.n_elem;
  arma::vec grad(p);

  #pragma omp parallel for default(none) \
      shared(M, block, row, col, symmetric, eps, grad, p) \
      firstprivate(logLik) \
      schedule(static)
  for (std::size_t k = 0; k < p; ++k) {
    LMSModel Mc_plus  = M.thread_clone();
    LMSModel Mc_minus = M.thread_clone();

    double& ti_plus  = lms_param(Mc_plus,  block[k], row[k], col[k]);
    double& ti_minus = lms_param(Mc_minus, block[k], row[k], col[k]);

    double* tj_plus  = nullptr;
    double* tj_minus = nullptr;

    if (symmetric[k] && row[k] != col[k]) {
      tj_plus  = &lms_param(Mc_plus,  block[k], col[k], row[k]);
      tj_minus = &lms_param(Mc_minus, block[k], col[k], row[k]);
    }

    ti_plus  += eps;
    ti_minus -= eps;
    if (tj_plus)  *tj_plus  += eps;
    if (tj_minus) *tj_minus -= eps;

    const double f_plus  = logLik(Mc_plus);
    const double f_minus = logLik(Mc_minus);

    grad[k] = (f_plus - f_minus) / (2.0 * eps);
  }

  return grad;
}

inline double completeLogLikFromModel(
    const LMSModel&  M,
    const arma::mat& V,
    const std::vector<arma::vec>& TGamma,
    const std::vector<std::vector<arma::vec>>& MeanPatterns,
    const std::vector<std::vector<arma::mat>>& CovPatterns,
    const std::vector<arma::uvec>& colidx,
    const arma::uvec n,
    const arma::uvec d,
    const int npatterns = 1) {

  const std::size_t J = V.n_rows;
  double ll = 0.0;

  for (std::size_t j = 0; j < J; j++) {

    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec& z  = V.row(j).t();   // view – no copy
    const arma::vec mu  = M.mu   (z);
    const arma::mat Sig = M.Sigma(z);

    for (int i = 0; i < npatterns; i++) {
      const arma::vec& nu = MeanPatterns[j][i];
      const arma::mat& S  = CovPatterns [j][i];
      const double tg = TGamma[j][i];

      if (tg <= DBL_MIN) continue;

      ll += totalDmvnWeighted(
        mu.elem(colidx[i]),
        Sig.submat(colidx[i], colidx[i]),
        nu, S, tg, d[i]);
    }
  }

  return ll;
}


// [[Rcpp::export]]
double completeLogLikLmsCpp(const Rcpp::List& modelR,
                            const Rcpp::List& P,
                            const Rcpp::List& quad,
                            const Rcpp::List& colidxR,
                            const arma::uvec& n,
                            const arma::uvec& d,
                            const int npatterns = 1) {
  const LMSModel model = LMSModel(modelR);
  const arma::mat Z = Rcpp::as<arma::mat>(quad["n"]).t(); // transpose so we can use column order vectors

  const arma::mat V       = Rcpp::as<arma::mat>(P["V"]);
  const auto      TGamma  = as_vec_of_vec(P["tgamma"]);
  const auto      Mean    = as_vec_of_vec_of_vec(P["mean"]);
  const auto      Cov     = as_vec_of_vec_of_mat(P["cov"]);
  const auto      colidx  = as_vec_of_uvec(colidxR);

  const Rcpp::List info   = modelR["info"];

  return completeLogLikFromModel(model, V, TGamma, Mean, Cov,
                                 colidx, n, d, npatterns);
}


// [[Rcpp::export]]
arma::vec gradLogLikLmsCpp(const Rcpp::List& modelR,
                           const Rcpp::List& P,
                           const arma::uvec& block,
                           const arma::uvec& row,
                           const arma::uvec& col,
                           const arma::uvec& symmetric,
                           const Rcpp::List& colidxR,
                           const arma::uvec& n,
                           const arma::uvec& d,
                           const int         npatterns = 1,
                           const double      eps = 1e-6,
                           const int         ncores = 1L) {
  LMSModel M(modelR);

  const arma::mat V       = Rcpp::as<arma::mat>(P["V"]);
  const auto      TGamma  = as_vec_of_vec(P["tgamma"]);
  const auto      Mean    = as_vec_of_vec_of_vec(P["mean"]);
  const auto      Cov     = as_vec_of_vec_of_mat(P["cov"]);
  const auto      colidx  = as_vec_of_uvec(colidxR);

  const arma::uvec idxAlpha      = arma::find(block == 8u);
  const arma::uvec idxBeta0      = arma::find(block == 9u);
  const arma::uvec idxTauX       = arma::find(block == 2u);
  const arma::uvec idxLambdaX    = arma::find(block == 0u);
  const arma::uvec idxA          = arma::find(block == 6u);
  const arma::uvec idxThetaDelta = arma::find(block == 4u);
  const arma::uvec idxPsi        = arma::find(block == 7u);
  const arma::uvec idxGammaXi    = arma::find(block == 10u);
  const arma::uvec idxGammaEta   = arma::find(block == 11u);
  const arma::uvec idxOmegaXiXi  = arma::find(block == 12u);
  const arma::uvec idxOmegaEtaXi = arma::find(block == 13u);
  arma::vec grad(block.n_elem, arma::fill::zeros);

  std::vector<arma::uword> fdPositions;
  fdPositions.reserve(block.n_elem);
  for (arma::uword pos = 0; pos < block.n_elem; ++pos) {
    const unsigned blk = block[pos];
    const bool analytic = (blk == 8u) || (blk == 9u) ||
                          (blk == 2u) || (blk == 0u) ||
                          (blk == 4u) || (blk == 7u) ||
                          (blk == 6u) || (blk == 10u) ||
                          (blk == 11u) || (blk == 12u) ||
                          (blk == 13u);
    if (!analytic) fdPositions.push_back(pos);
  }

  if (!fdPositions.empty()) {
    arma::uvec blockFD(fdPositions.size());
    arma::uvec rowFD(fdPositions.size());
    arma::uvec colFD(fdPositions.size());
    arma::uvec symFD(fdPositions.size());

    for (arma::uword i = 0; i < fdPositions.size(); ++i) {
      const arma::uword pos = fdPositions[i];
      blockFD[i] = block[pos];
      rowFD[i]   = row[pos];
      colFD[i]   = col[pos];
      symFD[i]   = symmetric[pos];
    }

    auto comp_ll = [&](LMSModel& mod) -> double {
      return completeLogLikFromModel(mod, V, TGamma, Mean, Cov,
                                     colidx, n, d, npatterns);
    };

    arma::vec gradFD = gradientCentralFD(M, comp_ll, blockFD, rowFD, colFD, symFD, eps, ncores);
    for (arma::uword i = 0; i < fdPositions.size(); ++i)
      grad[fdPositions[i]] = gradFD[i];
  }

  // Analytic gradient for alpha block
  if (!idxAlpha.is_empty()) {
    const arma::mat gradAlpha =
        gradAlphaComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxAlpha.n_elem; ++k) {
      const std::size_t pos = idxAlpha[k];
      grad[pos] = gradAlpha(row[pos], col[pos]);
    }
  }

  // Analytic gradient for beta0 block
  if (!idxBeta0.is_empty()) {
    const arma::mat gradBeta =
        gradBeta0Complete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxBeta0.n_elem; ++k) {
      const std::size_t pos = idxBeta0[k];
      grad[pos] = gradBeta(row[pos], col[pos]);
    }
  }

  // Analytic gradient for tauX block
  if (!idxTauX.is_empty()) {
    const arma::mat gradTau =
        gradTauXComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxTauX.n_elem; ++k) {
      const std::size_t pos = idxTauX[k];
      grad[pos] = gradTau(row[pos], col[pos]);
    }
  }

  // Analytic gradient for lambdaX block
  if (!idxLambdaX.is_empty()) {
    const arma::mat gradLambda =
        gradLambdaXComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxLambdaX.n_elem; ++k) {
      const std::size_t pos = idxLambdaX[k];
      grad[pos] = gradLambda(row[pos], col[pos]);
    }
  }

  // Analytic gradient for psi block
  if (!idxPsi.is_empty()) {
    const arma::mat gradPsi =
        gradPsiComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxPsi.n_elem; ++k) {
      const std::size_t pos = idxPsi[k];
      grad[pos] = gradPsi(row[pos], col[pos]);
    }
  }

  // Analytic gradient for thetaDelta block
  if (!idxThetaDelta.is_empty()) {
    const arma::mat gradTheta =
        gradThetaDeltaComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxThetaDelta.n_elem; ++k) {
      const std::size_t pos = idxThetaDelta[k];
      grad[pos] = gradTheta(row[pos], col[pos]);
    }
  }

  // Analytic gradient for A block
  if (!idxA.is_empty()) {
    const arma::mat gradA =
        gradAComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxA.n_elem; ++k) {
      const std::size_t pos = idxA[k];
      grad[pos] = gradA(row[pos], col[pos]);
    }
  }

  // Analytic gradient for gammaXi block
  if (!idxGammaXi.is_empty()) {
    const arma::mat gradGx =
        gradGammaXiComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxGammaXi.n_elem; ++k) {
      const std::size_t pos = idxGammaXi[k];
      grad[pos] = gradGx(row[pos], col[pos]);
    }
  }

  // Analytic gradient for gammaEta block
  if (!idxGammaEta.is_empty()) {
    const arma::mat gradGe =
        gradGammaEtaComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxGammaEta.n_elem; ++k) {
      const std::size_t pos = idxGammaEta[k];
      grad[pos] = gradGe(row[pos], col[pos]);
    }
  }

  // Analytic gradient for omegaXiXi block
  if (!idxOmegaXiXi.is_empty()) {
    const arma::mat gradOxx =
        gradOmegaXiXiComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxOmegaXiXi.n_elem; ++k) {
      const std::size_t pos = idxOmegaXiXi[k];
      grad[pos] = gradOxx(row[pos], col[pos]);
    }
  }

  // Analytic gradient for omegaEtaXi block
  if (!idxOmegaEtaXi.is_empty()) {
    const arma::mat gradOex =
        gradOmegaEtaXiComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxOmegaEtaXi.n_elem; ++k) {
      const std::size_t pos = idxOmegaEtaXi[k];
      grad[pos] = gradOex(row[pos], col[pos]);
    }
  }

  return grad;
}


// [[Rcpp::export]]
arma::vec gradLogLikLmsCpp_fd(const Rcpp::List& modelR,
                              const Rcpp::List& P,
                              const arma::uvec& block,
                              const arma::uvec& row,
                              const arma::uvec& col,
                              const arma::uvec& symmetric,
                              const Rcpp::List& colidxR,
                              const arma::uvec& n,
                              const arma::uvec& d,
                              const int         npatterns = 1,
                              const double      eps = 1e-6,
                              const int         ncores = 1L) {
  LMSModel M(modelR);

  const arma::mat V       = Rcpp::as<arma::mat>(P["V"]);
  const auto      TGamma  = as_vec_of_vec(P["tgamma"]);
  const auto      Mean    = as_vec_of_vec_of_vec(P["mean"]);
  const auto      Cov     = as_vec_of_vec_of_mat(P["cov"]);
  const auto      colidx  = as_vec_of_uvec(colidxR);

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(mod, V, TGamma, Mean, Cov,
                                   colidx, n, d, npatterns);
  };

  return gradientFD(M, comp_ll, block, row, col, symmetric, eps, ncores);
}


// [[Rcpp::export]]
arma::vec gradLogLikLmsCpp_cd(const Rcpp::List& modelR,
                              const Rcpp::List& P,
                              const arma::uvec& block,
                              const arma::uvec& row,
                              const arma::uvec& col,
                              const arma::uvec& symmetric,
                              const Rcpp::List& colidxR,
                              const arma::uvec& n,
                              const arma::uvec& d,
                              const int         npatterns = 1,
                              const double      eps = 1e-6,
                              const int         ncores = 1L) {
  LMSModel M(modelR);

  const arma::mat V       = Rcpp::as<arma::mat>(P["V"]);
  const auto      TGamma  = as_vec_of_vec(P["tgamma"]);
  const auto      Mean    = as_vec_of_vec_of_vec(P["mean"]);
  const auto      Cov     = as_vec_of_vec_of_mat(P["cov"]);
  const auto      colidx  = as_vec_of_uvec(colidxR);

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(mod, V, TGamma, Mean, Cov,
                                   colidx, n, d, npatterns);
  };

  return gradientCentralFD(M, comp_ll, block, row, col, symmetric, eps, ncores);
}


inline double observedLogLikFromModel(const LMSModel&  M,
                                      const arma::mat& V,
                                      const arma::vec& w,
                                      const arma::vec& samplingWeights,
                                      const std::vector<arma::mat>& data,
                                      const std::vector<arma::uvec>& colidx,
                                      const arma::uvec n,
                                      const int npatterns = 1,
                                      const int ncores = 1) {
  const std::size_t Q = V.n_rows;

  arma::vec density = arma::zeros<arma::vec>(arma::sum(n));

  for (std::size_t i = 0; i < Q; ++i) {
    if (w[i] <= DBL_MIN) continue;

    const arma::vec z   = V.row(i).t();
    const arma::vec mu  = M.mu   (z);
    const arma::mat Sig = M.Sigma(z);

    int offset = 0L;
    for (int j = 0; j < npatterns; j++) {
      const int end = offset + n[j] - 1L;

      density.subvec(offset, end) +=
        dmvnfast(data[j],
                 mu.elem(colidx[j]),
                 Sig.submat(colidx[j], colidx[j]),
                 false, ncores, false) * w[i];

      offset = end + 1L;
    }
  }

  return arma::sum(samplingWeights * arma::log(density));
}


// [[Rcpp::export]]
arma::vec gradObsLogLikLmsCpp(const Rcpp::List& modelR,
                              const Rcpp::List& dataR,
                              const Rcpp::List& colidxR,
                              const Rcpp::List& P,
                              const arma::uvec& block,
                              const arma::uvec& row,
                              const arma::uvec& col,
                              const arma::uvec& symmetric,
                              const arma::uvec& n,
                              const double      eps       = 1e-6,
                              const int         npatterns = 1L,
                              const int         ncores    = 1L) {
  LMSModel M(modelR);

  const arma::mat V = Rcpp::as<arma::mat>(P["V"]);
  const arma::vec w = Rcpp::as<arma::vec>(P["w"]);
  const arma::vec samplingWeights = Rcpp::as<arma::vec>(P["sampling.weights"]);

  const auto colidx = as_vec_of_uvec(colidxR);
  const auto data   = as_vec_of_mat(dataR);

  auto obs_ll = [&](LMSModel& mod) -> double {
    return observedLogLikFromModel(mod, V, w, samplingWeights, data, colidx, n, npatterns, 1L); // single-threaded
  };

  return gradientFD(M, obs_ll, block, row, col, symmetric, eps, ncores); // multi-thread here instead
}


// [[Rcpp::export]]
double observedLogLikLmsCpp(const Rcpp::List& modelR,
                            const Rcpp::List& dataR,
                            const Rcpp::List& colidxR,
                            const Rcpp::List& P,
                            const arma::uvec& n,
                            const int npatterns = 1L,
                            const int ncores = 1L) {
  const LMSModel M = LMSModel(modelR);

  const arma::mat V = Rcpp::as<arma::mat>(P["V"]);
  const arma::vec w = Rcpp::as<arma::vec>(P["w"]);
  const arma::vec samplingWeights = Rcpp::as<arma::vec>(P["sampling.weights"]);

  const auto colidx = as_vec_of_uvec(colidxR);
  const auto data   = as_vec_of_mat(dataR);

  return observedLogLikFromModel(M, V, w, samplingWeights, data, colidx, n, npatterns, ncores);
}


arma::mat gradAlphaComplete(const LMSModel&  M,
                            const arma::mat& V,
                            const std::vector<arma::vec>& TGamma,
                            const std::vector<std::vector<arma::vec>>& MeanPatterns,
                            const std::vector<std::vector<arma::mat>>& CovPatterns,
                            const std::vector<arma::uvec>& colidx,
                            const int npatterns) {
  arma::mat gradAlpha(M.a.n_rows, M.a.n_cols, arma::fill::zeros);

  const std::size_t J = V.n_rows;
  const unsigned numEta = M.Ie.n_rows;

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv = arma::inv(B);
    const arma::vec rhs = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;

    const arma::vec latent = arma::join_cols(muXi, muEta);
    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;

      const arma::uvec& idx = colidx[i];
      const arma::vec& nu = MeanPatterns[j][i];
      const arma::mat& S  = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);

      gradMu.elem(idx) += tg * linv_diff;
    }

    if (arma::all(gradMu == 0.0)) continue;

    arma::vec gradLatent = M.lX.t() * gradMu;
    arma::vec gradMuEta = gradLatent.subvec(M.numXis, M.numXis + numEta - 1);
    arma::vec gradAlphaNode = Binv.t() * gradMuEta;

    gradAlpha += arma::reshape(gradAlphaNode, gradAlpha.n_rows, gradAlpha.n_cols);
  }

  return gradAlpha;
}


arma::mat gradBeta0Complete(const LMSModel&  M,
                            const arma::mat& V,
                            const std::vector<arma::vec>& TGamma,
                            const std::vector<std::vector<arma::vec>>& MeanPatterns,
                            const std::vector<std::vector<arma::mat>>& CovPatterns,
                            const std::vector<arma::uvec>& colidx,
                            const int npatterns) {
  arma::mat gradBeta(M.beta0.n_rows, M.beta0.n_cols, arma::fill::zeros);

  const std::size_t J = V.n_rows;
  const unsigned numEta = M.Ie.n_rows;
  const arma::mat Oi = make_Oi(M.k, M.numXis);

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv = arma::inv(B);
    const arma::vec rhs = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;
    const arma::mat Eta = Binv * (M.Gx * M.A + kronZ.t() * M.Oxx * M.A);
    const arma::vec latent = arma::join_cols(muXi, muEta);
    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);
    std::vector<arma::mat> gradSigmaPatterns(npatterns);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) {
        gradSigmaPatterns[i].reset();
        continue;
      }

      const arma::uvec& idx = colidx[i];
      const arma::vec& nu = MeanPatterns[j][i];
      const arma::mat& S  = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) {
        gradSigmaPatterns[i].reset();
        continue;
      }

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);

      gradMu.elem(idx) += tg * linv_diff;

      arma::mat SigInv = arma::inv_sympd(SigSel);
      arma::mat diffOuter = diff * diff.t();
      arma::mat termS = SigInv * S * SigInv;
      arma::mat termDiff = SigInv * diffOuter * SigInv;
      gradSigmaPatterns[i] = 0.5 * (termS + tg * termDiff - tg * SigInv);
    }

    if (arma::all(gradMu == 0.0) && std::all_of(gradSigmaPatterns.begin(), gradSigmaPatterns.end(),
        [](const arma::mat& m){ return m.is_empty(); })) continue;

    const arma::vec gradLatent = M.lX.t() * gradMu;
    const arma::uword lenXi = muXi.n_elem;
    const arma::vec gradXi = gradLatent.subvec(0, lenXi - 1);
    const arma::vec gradEta = gradLatent.subvec(lenXi, lenXi + numEta - 1);

    const arma::mat Q = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat varXi = M.A * Oi * M.A.t();
    const arma::uword dimXi = varXi.n_rows;
    const arma::uword dimEta = Eta.n_rows;
    const arma::uword dimLatent = dimXi + dimEta;

    for (arma::uword paramIdx = 0; paramIdx < lenXi; ++paramIdx) {
      arma::vec d_muXi(lenXi, arma::fill::zeros);
      d_muXi[paramIdx] = 1.0;

      const arma::mat d_kronZ = arma::kron(M.Ie, d_muXi);
      const arma::mat dB = -d_kronZ.t() * M.Oex;
      const arma::mat dBinv = -Binv * dB * Binv;

      arma::vec d_rhs = M.Gx * d_muXi +
                        d_kronZ.t() * M.Oxx * muXi +
                        kronZ.t() * M.Oxx * d_muXi;
      arma::vec d_muEta = Binv * (d_rhs - dB * muEta);

      arma::mat dQ = d_kronZ.t() * M.Oxx * M.A;
      arma::mat dEta = dBinv * Q + Binv * dQ;

      arma::mat dVarEta = dEta * Oi * Eta.t() + Eta * Oi * dEta.t() +
                          dBinv * M.Psi * Binv.t() + Binv * M.Psi * dBinv.t();
      arma::mat dCovXiEta = M.A * Oi * dEta.t();

      arma::mat dVcov(dimLatent, dimLatent, arma::fill::zeros);
      dVcov.submat(0, dimXi, dimXi - 1, dimLatent - 1) = dCovXiEta;
      dVcov.submat(dimXi, 0, dimLatent - 1, dimXi - 1) = dCovXiEta.t();
      dVcov.submat(dimXi, dimXi, dimLatent - 1, dimLatent - 1) = dVarEta;

      arma::mat dSigma = M.lX * dVcov * M.lX.t();

      double cov_contrib = 0.0;
      for (int i = 0; i < npatterns; ++i) {
        if (gradSigmaPatterns[i].is_empty()) continue;
        const arma::uvec& idx = colidx[i];
        if (idx.is_empty()) continue;
        arma::mat dSigmaSub = dSigma.submat(idx, idx);
        cov_contrib += arma::accu(gradSigmaPatterns[i] % dSigmaSub);
      }

      const double mean_contrib = gradXi[paramIdx] + arma::dot(gradEta, d_muEta);
      const arma::uword r = paramIdx % M.beta0.n_rows;
      const arma::uword c = paramIdx / M.beta0.n_rows;
      gradBeta(r, c) += mean_contrib + cov_contrib;
    }
  }

  return gradBeta;
}


arma::mat gradTauXComplete(const LMSModel&  M,
                           const arma::mat& V,
                           const std::vector<arma::vec>& TGamma,
                           const std::vector<std::vector<arma::vec>>& MeanPatterns,
                           const std::vector<std::vector<arma::mat>>& CovPatterns,
                           const std::vector<arma::uvec>& colidx,
                           const int npatterns) {
  arma::mat gradTauX(M.tX.n_rows, M.tX.n_cols, arma::fill::zeros);

  const std::size_t J = V.n_rows;

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec mu = M.mu(z);
    const arma::mat Sigma = M.Sigma(z);

    arma::vec gradMu(mu.n_elem, arma::fill::zeros);

    for (int i = 0; i < npatterns; ++i) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;

      const arma::uvec& idx = colidx[i];
      const arma::vec& nu = MeanPatterns[j][i];
      const arma::mat& S  = CovPatterns [j][i];

      arma::mat SigSel = Sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, SigSel, "lower")) continue;

      const arma::vec diff = nu - mu.elem(idx);
      arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
      arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);

      gradMu.elem(idx) += tg * linv_diff;
    }

    if (arma::all(gradMu == 0.0)) continue;

    for (arma::uword c = 0; c < M.tX.n_cols; ++c) {
      for (arma::uword r = 0; r < M.tX.n_rows; ++r) {
        const arma::uword idx = r + c * M.tX.n_rows;
        if (idx < gradMu.n_elem)
          gradTauX(r, c) += gradMu(idx);
      }
    }
  }

  return gradTauX;
}




inline arma::vec getParams(const LMSModel& M,
                            const arma::uvec& block,
                            const arma::uvec& row,
                            const arma::uvec& col) {
  const std::size_t p = block.n_elem;
  arma::vec pars(p);
  for (std::size_t k = 0; k < p; ++k)
    pars[k] = lms_param(const_cast<LMSModel&>(M),
        block[k], row[k], col[k]);
  return pars;
}


inline void setParams(LMSModel&         M,
                       const arma::uvec& block,
                       const arma::uvec& row,
                       const arma::uvec& col,
                       const arma::uvec& symmetric,
                       const arma::vec&  vals) {
  const std::size_t p = block.n_elem;

  for (std::size_t k = 0; k < p; ++k) {
    double& ti = lms_param(M, block[k], row[k], col[k]);
    ti = vals[k];

    if (symmetric[k] && row[k] != col[k])
      lms_param(M, block[k], col[k], row[k]) = vals[k];
  }
}
template<class F>
Rcpp::List fdHessQuadraticFit(LMSModel&         M,
                                F&&               fun,
                                const arma::uvec& block,
                                const arma::uvec& row,
                                const arma::uvec& col,
                                const arma::uvec& symmetric,
                                const arma::vec&  base,
                                const arma::vec&  incr,
                                const int         ncores) {
  const std::size_t p = block.n_elem;

  // Build Koschal displacement matrix
  std::vector< arma::vec > disp;
  disp.reserve(1 + 2*p + (p*(p-1))/2);
  disp.emplace_back(arma::zeros<arma::vec>(p));
  for (std::size_t i = 0; i < p; ++i) {
    arma::vec v = arma::zeros<arma::vec>(p);
    v[i] =  1; disp.push_back(v);
    v[i] = -1; disp.push_back(v);
  }
  for (std::size_t i = 0; i < p-1; ++i)
    for (std::size_t j = i+1; j < p; ++j) {
      arma::vec v = arma::zeros<arma::vec>(p);
      v[i] = v[j] = 1;
      disp.push_back(v);
    }
  const std::size_t m = disp.size();

  // Evaluate fun at all design points (parallel)
  arma::vec y(m);
#pragma omp parallel for default(none) \
  shared(M, disp, m, block, row, col, symmetric, base, incr, y) \
  firstprivate(fun) schedule(static)
  for (std::size_t k = 0; k < m; ++k) {
    LMSModel Mc = M.thread_clone();
    setParams(Mc, block, row, col, symmetric, base + disp[k] % incr);
    y[k] = fun(Mc);
  }

  // Restore baseline
  setParams(M, block, row, col, symmetric, base);

  // Build design matrix
  const std::size_t q = 1 + 2*p + (p*(p-1))/2;
  arma::mat X(m, q, arma::fill::ones);
  std::size_t col_id = 1;
  for (std::size_t j = 0; j < p; ++j, ++col_id)
    for (std::size_t k = 0; k < m; ++k)
      X(k, col_id) = disp[k][j];
  for (std::size_t j = 0; j < p; ++j, ++col_id)
    for (std::size_t k = 0; k < m; ++k)
      X(k, col_id) = std::pow(disp[k][j], 2);
  for (std::size_t i = 0; i < p-1; ++i)
    for (std::size_t j = i+1; j < p; ++j, ++col_id)
      for (std::size_t k = 0; k < m; ++k)
        X(k, col_id) = disp[k][i] * disp[k][j];

  // frac scaling
  arma::vec frac(q, arma::fill::ones);
  for (std::size_t j = 0; j < p; ++j)              frac[1 + j]     = incr[j];
  for (std::size_t j = 0; j < p; ++j)              frac[1 + p + j] = incr[j]*incr[j];
  col_id = 1 + 2*p;
  for (std::size_t i = 0; i < p-1; ++i)
    for (std::size_t j = i+1; j < p; ++j, ++col_id)
      frac[col_id] = incr[i] * incr[j];

  arma::vec coef = arma::solve(X, y) / frac;

  arma::vec grad = coef.subvec(1, p);
  arma::mat Hess(p, p, arma::fill::zeros);
  for (std::size_t j = 0; j < p; ++j)
    Hess(j, j) = 2.0 * coef[1 + p + j];
  col_id = 1 + 2*p;
  for (std::size_t i = 0; i < p-1; ++i)
    for (std::size_t j = i+1; j < p; ++j, ++col_id) {
      Hess(i,j) = coef[col_id];
      Hess(j,i) = coef[col_id];
    }

  return Rcpp::List::create(
    Rcpp::Named("mean")     = coef[0],
    Rcpp::Named("gradient") = grad,
    Rcpp::Named("Hessian")  = Hess
  );
}


template<class F>
Rcpp::List fdHessFullFd(LMSModel&         M,
                        F&&               fun,
                        const arma::uvec& block,
                        const arma::uvec& row,
                        const arma::uvec& col,
                        const arma::uvec& symmetric,
                        const arma::vec&  base,
                        const arma::vec&  incr,
                        const int         ncores) {
  const std::size_t p = block.n_elem;
  const std::size_t npairs = (p>1) ? (p*(p-1))/2 : 0;
  const std::size_t m = 1 + 2*p + 4*npairs;

  // Index helper for pairs
  auto pairIndex = [p](std::size_t i, std::size_t j) -> std::size_t {
    return (i*(2*p - i - 1))/2 + (j - i - 1);
  };

  // Build displacements
  std::vector< arma::vec > disp;
  disp.reserve(m);
  disp.emplace_back(arma::zeros<arma::vec>(p)); // origin
  const std::size_t idx0 = 0;

  std::vector<std::size_t> idx_ip(p), idx_im(p);
  for (std::size_t i=0; i<p; ++i) {
    arma::vec v = arma::zeros<arma::vec>(p);
    v[i]= 1; idx_ip[i]=disp.size(); disp.push_back(v);
    v[i]=-1; idx_im[i]=disp.size(); disp.push_back(v);
  }

  std::vector<std::size_t> idx_pp(npairs), idx_pm(npairs),
                           idx_mp(npairs), idx_mm(npairs);
  if (p>1) {
    for (std::size_t i=0; i<p-1; ++i)
      for (std::size_t j=i+1; j<p; ++j) {
        std::size_t k = pairIndex(i,j);
        arma::vec v = arma::zeros<arma::vec>(p);
        v[i]= 1; v[j]= 1; idx_pp[k]=disp.size(); disp.push_back(v);
        v[i]= 1; v[j]=-1; idx_pm[k]=disp.size(); disp.push_back(v);
        v[i]=-1; v[j]= 1; idx_mp[k]=disp.size(); disp.push_back(v);
        v[i]=-1; v[j]=-1; idx_mm[k]=disp.size(); disp.push_back(v);
      }
  }

  // Evaluate fun (parallel)
  arma::vec y(disp.size());
#pragma omp parallel for default(none) \
  shared(M, disp, block, row, col, symmetric, base, incr, y) \
  firstprivate(fun) schedule(static)
  for (std::size_t k=0; k<disp.size(); ++k) {
    LMSModel Mc = M.thread_clone();
    setParams(Mc, block, row, col, symmetric, base + disp[k] % incr);
    y[k] = fun(Mc);
  }
  setParams(M, block, row, col, symmetric, base);

  // Assemble gradient/Hessian
  arma::vec grad(p, arma::fill::zeros);
  arma::mat Hess(p, p, arma::fill::zeros);
  const double f0 = y[idx0];

  for (std::size_t i=0; i<p; ++i) {
    double hi = incr[i];
    double f_ip = y[idx_ip[i]];
    double f_im = y[idx_im[i]];
    grad[i]  = (f_ip - f_im) / (2.0*hi);
    Hess(i,i)= (f_ip + f_im - 2.0*f0) / (hi*hi);
  }

  if (p>1) {
    for (std::size_t i=0; i<p-1; ++i) {
      double hi = incr[i];
      for (std::size_t j=i+1; j<p; ++j) {
        double hj = incr[j];
        std::size_t k = pairIndex(i,j);
        double fpp=y[idx_pp[k]], fpm=y[idx_pm[k]],
               fmp=y[idx_mp[k]], fmm=y[idx_mm[k]];
        double hij = (fpp - fpm - fmp + fmm) / (4.0*hi*hj);
        Hess(i,j)=hij; Hess(j,i)=hij;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mean")     = f0,
    Rcpp::Named("gradient") = grad,
    Rcpp::Named("Hessian")  = Hess
  );
}


// ======================================================
// Dispatcher
// ======================================================
template<class F>
Rcpp::List fdHessCpp(LMSModel&         M,
                     F&&               fun,
                     const arma::uvec& block,
                     const arma::uvec& row,
                     const arma::uvec& col,
                     const arma::uvec& symmetric,
                     const double      relStep   = 1e-6,
                     const double      minAbsPar = 0.0,
                     const int         ncores    = 1L) {
  ThreadSetter ts(ncores);

  const std::size_t p = block.n_elem;
  const arma::vec base = getParams(M, block, row, col);
  const arma::vec incr =
      arma::max(arma::abs(base), arma::vec(p).fill(minAbsPar)) * relStep;

  // Switching heuristics
  constexpr std::size_t P_SWITCH        = 120;
  constexpr std::size_t MEM_LIMIT_BYTES = 3ull << 30;

  auto m_ls = 1 + 2*p + (p*(p-1))/2;
  auto bytes_X = (unsigned long long)m_ls * (unsigned long long)m_ls *
                 (unsigned long long)sizeof(double);

  bool useFullFd = (p >= P_SWITCH) || (bytes_X > MEM_LIMIT_BYTES);

  if (!useFullFd)
    return fdHessQuadraticFit(M, std::forward<F>(fun), block, row, col,
        symmetric, base, incr, ncores);
  else
    return fdHessFullFd(M, std::forward<F>(fun), block, row, col,
        symmetric, base, incr, ncores);
}


// [[Rcpp::export]]
Rcpp::List hessObsLogLikLmsCpp(const Rcpp::List& modelR,
                               const Rcpp::List& dataR,
                               const Rcpp::List& P,
                               const arma::uvec& block,
                               const arma::uvec& row,
                               const arma::uvec& col,
                               const arma::uvec& symmetric,
                               const Rcpp::List& colidxR,
                               const arma::uvec& n,
                               const int         npatterns = 1L,
                               const double      relStep = 1e-6,
                               const double      minAbs  = 0.0,
                               const int         ncores  = 1L) {
    LMSModel M(modelR);

    const arma::mat V = Rcpp::as<arma::mat>(P["V"]);
    const arma::vec w = Rcpp::as<arma::vec>(P["w"]);
    const arma::vec samplingWeights = Rcpp::as<arma::vec>(P["sampling.weights"]);

    const auto colidx = as_vec_of_uvec(colidxR);
    const auto data   = as_vec_of_mat(dataR);

    auto obs_ll = [&](LMSModel& mod) -> double {
        return observedLogLikFromModel(mod, V, w, samplingWeights, data, colidx,
                                       n, npatterns, 1L); // single-threaded
    };

    return fdHessCpp(M, obs_ll, block, row, col, symmetric,
        relStep, minAbs, ncores); // multi-threaded
}


// [[Rcpp::export]]
Rcpp::List hessCompLogLikLmsCpp(const Rcpp::List& modelR,
                                const Rcpp::List& P,
                                const arma::uvec& block,
                                const arma::uvec& row,
                                const arma::uvec& col,
                                const arma::uvec& symmetric,
                                const Rcpp::List& colidxR,
                                const arma::uvec& n,
                                const arma::uvec& d,
                                const int         npatterns = 1,
                                const double      relStep   = 1e-6,
                                const double      minAbs    = 0.0,
                                const int         ncores    = 1L) {
  LMSModel M(modelR);

  const arma::mat  V       = Rcpp::as<arma::mat>(P["V"]);
  const auto       TGamma  = as_vec_of_vec(P["tgamma"]);
  const auto       Mean    = as_vec_of_vec_of_vec(P["mean"]);
  const auto       Cov     = as_vec_of_vec_of_mat(P["cov"]);
  const auto       colidx  = as_vec_of_uvec(colidxR);

  const Rcpp::List info   = modelR["info"];

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(mod, V, TGamma, Mean, Cov,
                                   colidx, n, d, npatterns); // single-threaded
  };

  return fdHessCpp(M, comp_ll, block, row, col, symmetric,
      relStep, minAbs, ncores); // multi-threaded
}
