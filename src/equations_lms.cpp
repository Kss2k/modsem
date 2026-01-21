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


struct AlphaDerivatives {
  arma::mat grad;
  arma::mat hess;
};

struct TauXDerivatives {
  arma::mat grad;
  arma::mat hess;
};

struct Beta0Derivatives {
  arma::mat grad;
  arma::mat hess;
};

struct LambdaXDerivatives {
  arma::mat grad;
  arma::mat hess;
};

struct ThetaDeltaDerivatives {
  arma::mat grad;
  arma::mat hess;
};

struct PsiDerivatives {
  arma::mat grad;
  arma::mat hess;
};

struct GammaXiDerivatives {
  arma::mat grad;
  arma::mat hess;
};

struct GammaEtaDerivatives {
  arma::mat grad;
  arma::mat hess;
};


struct CompleteCaseDerivatives {
  arma::vec gradMu;
  arma::mat curvMu;
  arma::mat gradSigma;
  std::vector<arma::mat> gradSigmaPatterns;
  std::vector<arma::mat> sigmaInvPatterns;
  std::vector<arma::vec> diffPatterns;
  std::vector<arma::vec> precDiffPatterns;
  std::vector<double> patternWeights;

  bool hasMeanGrad() const {
    return !gradMu.is_empty() && !gradMu.is_zero();
  }

  bool hasMeanCurv() const {
    return !curvMu.is_empty() && !curvMu.is_zero();
  }

  bool hasCovGrad() const {
    return !gradSigma.is_empty() && !gradSigma.is_zero();
  }
};


CompleteCaseDerivatives computeCompleteCaseDerivatives(
    const arma::vec& mu,
    const arma::mat& Sigma,
    const arma::vec& patternWeights,
    const std::vector<arma::vec>& meanPatterns,
    const std::vector<arma::mat>& covPatterns,
    const std::vector<arma::uvec>& colidx,
    const int npatterns) {

  CompleteCaseDerivatives out;
  out.gradMu = arma::vec(mu.n_elem, arma::fill::zeros);
  out.curvMu = arma::mat(mu.n_elem, mu.n_elem, arma::fill::zeros);
  out.gradSigma = arma::mat(Sigma.n_rows, Sigma.n_cols, arma::fill::zeros);
  out.gradSigmaPatterns.resize(npatterns);
  out.sigmaInvPatterns.resize(npatterns);
  out.diffPatterns.resize(npatterns);
  out.precDiffPatterns.resize(npatterns);
  out.patternWeights.resize(npatterns, 0.0);

  for (int i = 0; i < npatterns; ++i) {
    const bool hasWeight = patternWeights.n_elem > static_cast<arma::uword>(i);
    const double tg = hasWeight ? patternWeights[i] : 0.0;
    out.patternWeights[i] = tg;
    if (tg <= DBL_MIN) {
      out.gradSigmaPatterns[i].reset();
      out.sigmaInvPatterns[i].reset();
      out.diffPatterns[i].reset();
      out.precDiffPatterns[i].reset();
      continue;
    }

    if (static_cast<std::size_t>(i) >= colidx.size() ||
        static_cast<std::size_t>(i) >= meanPatterns.size() ||
        static_cast<std::size_t>(i) >= covPatterns.size()) {
      out.gradSigmaPatterns[i].reset();
      out.sigmaInvPatterns[i].reset();
      out.diffPatterns[i].reset();
      out.precDiffPatterns[i].reset();
      continue;
    }

    const arma::uvec& idx = colidx[i];
    if (idx.is_empty()) {
      out.gradSigmaPatterns[i].reset();
      out.sigmaInvPatterns[i].reset();
      out.diffPatterns[i].reset();
      out.precDiffPatterns[i].reset();
      continue;
    }

    const arma::vec& nu = meanPatterns[i];
    const arma::mat& S  = covPatterns [i];
    if (nu.is_empty() || S.is_empty()) {
      out.gradSigmaPatterns[i].reset();
      out.sigmaInvPatterns[i].reset();
      out.diffPatterns[i].reset();
      out.precDiffPatterns[i].reset();
      continue;
    }

    if (nu.n_elem != idx.n_elem ||
        S.n_rows != idx.n_elem ||
        S.n_cols != idx.n_elem) {
      out.gradSigmaPatterns[i].reset();
      out.sigmaInvPatterns[i].reset();
      out.diffPatterns[i].reset();
      out.precDiffPatterns[i].reset();
      continue;
    }

    arma::mat SigSel = Sigma.submat(idx, idx);
    arma::mat L;
    if (!arma::chol(L, SigSel, "lower")) {
      out.gradSigmaPatterns[i].reset();
      out.sigmaInvPatterns[i].reset();
      out.diffPatterns[i].reset();
      out.precDiffPatterns[i].reset();
      continue;
    }

    const arma::vec diff = nu - mu.elem(idx);
    arma::vec y = arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast);
    arma::vec linv_diff = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);
    out.gradMu.elem(idx) += tg * linv_diff;

    arma::mat SigInv = arma::inv_sympd(SigSel);
    out.curvMu.submat(idx, idx) += tg * SigInv;

    arma::mat diffOuter = diff * diff.t();
    arma::mat termS = SigInv * S * SigInv;
    arma::mat termDiff = SigInv * diffOuter * SigInv;
    out.gradSigmaPatterns[i] = 0.5 * (termS + tg * termDiff - tg * SigInv);
    out.gradSigma.submat(idx, idx) += out.gradSigmaPatterns[i];
    out.sigmaInvPatterns[i] = SigInv;
    out.diffPatterns[i] = diff;
    out.precDiffPatterns[i] = linv_diff;
  }

  return out;
}

AlphaDerivatives alphaGradHessComplete(const LMSModel&  M,
                                       const arma::mat& V,
                                       const std::vector<arma::vec>& TGamma,
                                       const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                       const std::vector<std::vector<arma::mat>>& CovPatterns,
                                       const std::vector<arma::uvec>& colidx,
                                       const int npatterns);

TauXDerivatives tauXGradHessComplete(const LMSModel&  M,
                                     const arma::mat& V,
                                     const std::vector<arma::vec>& TGamma,
                                     const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                     const std::vector<std::vector<arma::mat>>& CovPatterns,
                                     const std::vector<arma::uvec>& colidx,
                                     const int npatterns);

arma::mat gradAlphaComplete(const LMSModel&  M,
                            const arma::mat& V,
                            const std::vector<arma::vec>& TGamma,
                            const std::vector<std::vector<arma::vec>>& MeanPatterns,
                            const std::vector<std::vector<arma::mat>>& CovPatterns,
                            const std::vector<arma::uvec>& colidx,
                            const int npatterns);


TauXDerivatives tauXGradHessComplete(const LMSModel&  M,
                                     const arma::mat& V,
                                     const std::vector<arma::vec>& TGamma,
                                     const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                     const std::vector<std::vector<arma::mat>>& CovPatterns,
                                     const std::vector<arma::uvec>& colidx,
                                     const int npatterns) {
  TauXDerivatives out;
  out.grad = arma::mat(M.tX.n_rows, M.tX.n_cols, arma::fill::zeros);
  out.hess = arma::mat(M.tX.n_elem, M.tX.n_elem, arma::fill::zeros);

  const std::size_t J = V.n_rows;
  const arma::uword nrows = M.tX.n_rows;
  const arma::uword ncols = M.tX.n_cols;
  const arma::uword nParams = M.tX.n_elem;

  auto mu_index = [&](arma::uword r, arma::uword c) -> arma::uword {
    return r + c * nrows;
  };

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const arma::vec mu = M.mu(z);
    const arma::mat Sigma = M.Sigma(z);

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const arma::mat& curvMu = cache.curvMu;

    const bool hasGrad = cache.hasMeanGrad();
    const bool hasCurv = cache.hasMeanCurv();
    if (!hasGrad && !hasCurv) continue;

    if (hasGrad) {
      for (arma::uword c = 0; c < ncols; ++c) {
        for (arma::uword r = 0; r < nrows; ++r) {
          const arma::uword idxMu = mu_index(r, c);
          if (idxMu >= gradMu.n_elem) continue;
          out.grad(r, c) += gradMu[idxMu];
        }
      }
    }

    if (hasCurv) {
      for (arma::uword p = 0; p < nParams; ++p) {
        const arma::uword rp = p % nrows;
        const arma::uword cp = p / nrows;
        const arma::uword idxP = mu_index(rp, cp);
        if (idxP >= curvMu.n_rows) continue;

        for (arma::uword q = 0; q < nParams; ++q) {
          const arma::uword rq = q % nrows;
          const arma::uword cq = q / nrows;
          const arma::uword idxQ = mu_index(rq, cq);
          if (idxQ >= curvMu.n_cols) continue;

          out.hess(p, q) -= curvMu(idxP, idxQ);
        }
      }
    }
  }

  return out;
}


arma::mat gradBeta0Complete(const LMSModel&  M,
                            const arma::mat& V,
                            const std::vector<arma::vec>& TGamma,
                            const std::vector<std::vector<arma::vec>>& MeanPatterns,
                            const std::vector<std::vector<arma::mat>>& CovPatterns,
                            const std::vector<arma::uvec>& colidx,
                            const int npatterns);

Beta0Derivatives beta0GradHessComplete(const LMSModel&  M,
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

LambdaXDerivatives lambdaXGradHessComplete(const LMSModel&  M,
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

struct ADerivatives {
  arma::mat grad;
  arma::mat hess;
};


arma::mat gradThetaDeltaComplete(const LMSModel&  M,
                                 const arma::mat& V,
                                 const std::vector<arma::vec>& TGamma,
                                 const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                 const std::vector<std::vector<arma::mat>>& CovPatterns,
                                 const std::vector<arma::uvec>& colidx,
                                 const int npatterns);

PsiDerivatives psiGradHessComplete(const LMSModel&  M,
                                   const arma::mat& V,
                                   const std::vector<arma::vec>& TGamma,
                                   const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                   const std::vector<std::vector<arma::mat>>& CovPatterns,
                                   const std::vector<arma::uvec>& colidx,
                                   const int npatterns);

ThetaDeltaDerivatives thetaDeltaGradHessComplete(const LMSModel&  M,
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

GammaXiDerivatives gammaXiGradHessComplete(const LMSModel&  M,
                                           const arma::mat& V,
                                           const std::vector<arma::vec>& TGamma,
                                           const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                           const std::vector<std::vector<arma::mat>>& CovPatterns,
                                           const std::vector<arma::uvec>& colidx,
                                           const int npatterns);

GammaEtaDerivatives gammaEtaGradHessComplete(const LMSModel&  M,
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

LambdaXDerivatives lambdaXGradHessComplete(const LMSModel&  M,
                                           const arma::mat& V,
                                           const std::vector<arma::vec>& TGamma,
                                           const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                           const std::vector<std::vector<arma::mat>>& CovPatterns,
                                           const std::vector<arma::uvec>& colidx,
                                           const int npatterns) {
  LambdaXDerivatives out;
  out.grad = arma::mat(M.lX.n_rows, M.lX.n_cols, arma::fill::zeros);
  out.hess = arma::mat(M.lX.n_elem, M.lX.n_elem, arma::fill::zeros);

  const arma::uword nrows = M.lX.n_rows;
  const arma::uword ncols = M.lX.n_cols;
  const arma::uword nParams = M.lX.n_elem;

  std::vector<arma::uword> paramRow(nParams);
  std::vector<arma::uword> paramCol(nParams);
  for (arma::uword c = 0; c < ncols; ++c) {
    for (arma::uword r = 0; r < nrows; ++r) {
      const arma::uword idx = r + c * nrows;
      paramRow[idx] = r;
      paramCol[idx] = c;
    }
  }

  std::vector<std::vector<int>> patternPos(npatterns);
  for (int i = 0; i < npatterns; ++i) {
    patternPos[i].assign(nrows, -1);
    const arma::uvec& idx = colidx[i];
    for (arma::uword pos = 0; pos < idx.n_elem; ++pos) {
      const arma::uword obs = idx[pos];
      if (obs < nrows) patternPos[i][obs] = static_cast<int>(pos);
    }
  }

  struct LambdaPatternDeriv {
    arma::vec dMu;
    arma::mat dSigma;
    int       muPos = -1;
  };

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
    const arma::mat LXvcov = M.lX * vcov;

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const std::vector<arma::mat>& gradSigmaPatterns = cache.gradSigmaPatterns;

    std::vector<std::vector<LambdaPatternDeriv>> paramPatterns(
        nParams, std::vector<LambdaPatternDeriv>(npatterns));

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::uword r = paramRow[p];
      const arma::uword c = paramCol[p];
      const double latentVal = latent[c];

      double mean_contrib = gradMu[r] * latentVal;
      double cov_contrib = 0.0;

      for (int i = 0; i < npatterns; ++i) {
        if (cache.patternWeights[i] <= DBL_MIN) continue;
        if (cache.sigmaInvPatterns[i].is_empty()) continue;
        const arma::uvec& idx = colidx[i];
        if (idx.is_empty()) continue;

        LambdaPatternDeriv deriv;
        const int pos = (r < nrows) ? patternPos[i][r] : -1;
        if (pos >= 0) {
          deriv.muPos = pos;
          deriv.dMu = arma::vec(idx.n_elem, arma::fill::zeros);
          deriv.dMu[pos] = latentVal;

          arma::vec w(idx.n_elem, arma::fill::zeros);
          for (arma::uword k = 0; k < idx.n_elem; ++k) {
            const arma::uword obs = idx[k];
            if (obs < LXvcov.n_rows) w[k] = LXvcov(obs, c);
          }

          deriv.dSigma = arma::mat(idx.n_elem, idx.n_elem, arma::fill::zeros);
          deriv.dSigma.row(pos) += w.t();
          deriv.dSigma.col(pos) += w;

          if (!gradSigmaPatterns[i].is_empty())
            cov_contrib += arma::accu(gradSigmaPatterns[i] % deriv.dSigma);
        }

        paramPatterns[p][i] = std::move(deriv);
      }

      out.grad(r, c) += mean_contrib + cov_contrib;
    }

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::uword r_p = paramRow[p];
      const arma::uword c_p = paramCol[p];
      const double latent_p = latent[c_p];

      for (arma::uword q = p; q < nParams; ++q) {
        const arma::uword r_q = paramRow[q];
        const arma::uword c_q = paramCol[q];

        double hessContrib = 0.0;

        for (int i = 0; i < npatterns; ++i) {
          if (cache.patternWeights[i] <= DBL_MIN) continue;
          if (cache.sigmaInvPatterns[i].is_empty()) continue;
          const arma::uvec& idx = colidx[i];
          if (idx.is_empty()) continue;

          const LambdaPatternDeriv& derivP = paramPatterns[p][i];
          const LambdaPatternDeriv& derivQ = paramPatterns[q][i];
          if (derivP.muPos < 0 && derivP.dSigma.is_empty()) continue;
          if (derivQ.muPos < 0 && derivQ.dSigma.is_empty()) continue;

          const arma::mat& A = cache.sigmaInvPatterns[i];
          const arma::vec& diff = cache.diffPatterns[i];
          const arma::mat& S = CovPatterns[j][i];
          const arma::mat& gradSigma = gradSigmaPatterns[i];
          const double tg = cache.patternWeights[i];

          arma::mat Cq = derivQ.dSigma.is_empty()
              ? arma::mat(A.n_rows, A.n_cols, arma::fill::zeros)
              : derivQ.dSigma;
          arma::vec Bq = derivQ.dMu.is_empty()
              ? arma::vec(A.n_rows, arma::fill::zeros)
              : derivQ.dMu;

          arma::mat dA = Cq.is_empty() ? arma::mat(A.n_rows, A.n_cols, arma::fill::zeros)
                                       : -A * Cq * A;
          arma::vec dDiff = Bq.is_empty() ? arma::vec(A.n_rows, arma::fill::zeros)
                                          : -Bq;

          double meanTerm = 0.0;
          if (derivP.muPos >= 0) {
            arma::vec dA_diff = dA * diff;
            arma::vec ABq = Bq.is_empty() ? arma::vec(A.n_rows, arma::fill::zeros)
                                          : A * Bq;
            meanTerm = tg * latent_p *
                ((derivP.muPos < static_cast<int>(dA_diff.n_elem) ? dA_diff[derivP.muPos] : 0.0) -
                 (derivP.muPos < static_cast<int>(ABq.n_elem) ? ABq[derivP.muPos] : 0.0));
          }

          arma::mat dGradSigma(A.n_rows, A.n_cols, arma::fill::zeros);
          if (!Cq.is_empty() || !Bq.is_empty()) {
            dGradSigma =
                0.5 * (dA * S * A + A * S * dA +
                       tg * (dA * diff * diff.t() * A +
                             A * diff * diff.t() * dA +
                             A * dDiff * diff.t() * A +
                             A * diff * dDiff.t() * A) -
                       tg * dA);
          }

          double covTerm = 0.0;
          if (!derivP.dSigma.is_empty())
            covTerm += arma::accu(dGradSigma % derivP.dSigma);

          if (!gradSigma.is_empty()) {
            arma::mat Cpq(A.n_rows, A.n_cols, arma::fill::zeros);
            if (derivP.muPos >= 0 && derivQ.muPos >= 0) {
              const double val = vcov(c_p, c_q);
              Cpq(derivP.muPos, derivQ.muPos) += val;
              Cpq(derivQ.muPos, derivP.muPos) += val;
            }
            covTerm += arma::accu(gradSigma % Cpq);
          }

          hessContrib += meanTerm + covTerm;
        }

        out.hess(p, q) += hessContrib;
        if (p != q) out.hess(q, p) += hessContrib;
      }
    }
  }

  return out;
}

arma::mat gradLambdaXComplete(const LMSModel&  M,
                              const arma::mat& V,
                              const std::vector<arma::vec>& TGamma,
                              const std::vector<std::vector<arma::vec>>& MeanPatterns,
                              const std::vector<std::vector<arma::mat>>& CovPatterns,
                              const std::vector<arma::uvec>& colidx,
                              const int npatterns) {
  return lambdaXGradHessComplete(M, V, TGamma, MeanPatterns, CovPatterns,
                                 colidx, npatterns).grad;
}


ThetaDeltaDerivatives thetaDeltaGradHessComplete(const LMSModel&  M,
                                                 const arma::mat& V,
                                                 const std::vector<arma::vec>& TGamma,
                                                 const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                                 const std::vector<std::vector<arma::mat>>& CovPatterns,
                                                 const std::vector<arma::uvec>& colidx,
                                                 const int npatterns) {
  ThetaDeltaDerivatives out;
  out.grad = arma::mat(M.d.n_rows, M.d.n_cols, arma::fill::zeros);
  out.hess = arma::mat(M.d.n_elem, M.d.n_elem, arma::fill::zeros);

  const arma::uword nrows = M.d.n_rows;
  const arma::uword ncols = M.d.n_cols;
  const arma::uword nParams = M.d.n_elem;

  std::vector<arma::uword> paramRow(nParams);
  std::vector<arma::uword> paramCol(nParams);
  for (arma::uword c = 0; c < ncols; ++c) {
    for (arma::uword r = 0; r < nrows; ++r) {
      const arma::uword idx = r + c * nrows;
      paramRow[idx] = r;
      paramCol[idx] = c;
    }
  }

  std::vector<std::vector<int>> patternPos(npatterns);
  for (int i = 0; i < npatterns; ++i) {
    patternPos[i].assign(nrows, -1);
    const arma::uvec& idx = colidx[i];
    for (arma::uword pos = 0; pos < idx.n_elem; ++pos) {
      const arma::uword obs = idx[pos];
      if (obs < nrows) patternPos[i][obs] = static_cast<int>(pos);
    }
  }

  struct ThetaPatternDeriv {
    arma::mat dSigma;
    bool active = false;
  };

  std::vector<std::vector<ThetaPatternDeriv>> paramPatterns(
      nParams, std::vector<ThetaPatternDeriv>(npatterns));

  for (arma::uword p = 0; p < nParams; ++p) {
    const arma::uword r = paramRow[p];
    const arma::uword c = paramCol[p];
    for (int i = 0; i < npatterns; ++i) {
      const arma::uvec& idx = colidx[i];
      if (idx.is_empty()) continue;
      const int pos_r = patternPos[i][r];
      const int pos_c = patternPos[i][c];
      if (pos_r < 0 || pos_c < 0) continue;

      ThetaPatternDeriv deriv;
      deriv.dSigma = arma::mat(idx.n_elem, idx.n_elem, arma::fill::zeros);
      deriv.dSigma(pos_r, pos_c) = 1.0;
      deriv.active = true;
      paramPatterns[p][i] = std::move(deriv);
    }
  }

  const std::size_t J = V.n_rows;

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z  = V.row(j).t();
    const arma::vec mu = M.mu(z);
    const arma::mat Sigma = M.Sigma(z);

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    if (!cache.hasCovGrad()) continue;

    const arma::uword rowLimit = std::min<arma::uword>(out.grad.n_rows, cache.gradSigma.n_rows);
    const arma::uword colLimit = std::min<arma::uword>(out.grad.n_cols, cache.gradSigma.n_cols);
    out.grad.submat(0, 0, rowLimit - 1, colLimit - 1) +=
        cache.gradSigma.submat(0, 0, rowLimit - 1, colLimit - 1);

    for (arma::uword p = 0; p < nParams; ++p) {
      for (arma::uword q = p; q < nParams; ++q) {
        double hessContrib = 0.0;

        for (int i = 0; i < npatterns; ++i) {
          if (cache.patternWeights[i] <= DBL_MIN) continue;
          if (cache.sigmaInvPatterns[i].is_empty()) continue;
          const arma::uvec& idx = colidx[i];
          if (idx.is_empty()) continue;

          const ThetaPatternDeriv& derivP = paramPatterns[p][i];
          const ThetaPatternDeriv& derivQ = paramPatterns[q][i];
          if (!derivP.active || !derivQ.active) continue;

          const arma::mat& A = cache.sigmaInvPatterns[i];
          const arma::vec& diff = cache.diffPatterns[i];
          const arma::mat& S = CovPatterns[j][i];
          const double tg = cache.patternWeights[i];

          arma::mat Cq = derivQ.dSigma;
          arma::mat dA = -A * Cq * A;

          arma::mat dGradSigma =
              0.5 * (dA * S * A + A * S * dA +
                     tg * (dA * diff * diff.t() * A +
                           A * diff * diff.t() * dA) -
                     tg * dA);

          hessContrib += arma::accu(dGradSigma % derivP.dSigma);
        }

        out.hess(p, q) += hessContrib;
        if (p != q) out.hess(q, p) += hessContrib;
      }
    }
  }

  return out;
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

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const arma::mat& gradSigma = cache.gradSigma;
    if (!cache.hasMeanGrad() && !cache.hasCovGrad()) continue;

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
        d_muXi[r] = (c < zVec.n_elem) ? zVec[c] : 0.0;
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

PsiDerivatives psiGradHessComplete(const LMSModel&  M,
                                   const arma::mat& V,
                                   const std::vector<arma::vec>& TGamma,
                                   const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                   const std::vector<std::vector<arma::mat>>& CovPatterns,
                                   const std::vector<arma::uvec>& colidx,
                                   const int npatterns) {
  PsiDerivatives out;
  out.grad = arma::mat(M.Psi.n_rows, M.Psi.n_cols, arma::fill::zeros);
  out.hess = arma::mat(M.Psi.n_elem, M.Psi.n_elem, arma::fill::zeros);

  const arma::uword dimEta = M.Psi.n_rows;
  if (dimEta == 0) return out;

  const arma::uword dimXi  = M.numXis;
  const arma::uword etaStart = dimXi;
  const arma::uword etaEnd   = etaStart + dimEta - 1;
  const arma::mat L_eta = M.lX.cols(etaStart, etaEnd);

  const arma::uword nParams = M.Psi.n_elem;
  std::vector<arma::uword> paramRow(nParams);
  std::vector<arma::uword> paramCol(nParams);
  for (arma::uword c = 0; c < M.Psi.n_cols; ++c) {
    for (arma::uword r = 0; r < M.Psi.n_rows; ++r) {
      const arma::uword idx = r + c * M.Psi.n_rows;
      paramRow[idx] = r;
      paramCol[idx] = c;
    }
  }

  const std::size_t J = V.n_rows;
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

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    if (!cache.hasCovGrad()) continue;

    const arma::mat& gradSigma = cache.gradSigma;
    const arma::mat Gvcov = M.lX.t() * gradSigma * M.lX;
    if (Gvcov.n_rows <= etaEnd || Gvcov.n_cols <= etaEnd) continue;
    const arma::mat GvarEta = Gvcov.submat(etaStart, etaStart, etaEnd, etaEnd);
    out.grad += Binv.t() * GvarEta * Binv;

    std::vector<std::vector<arma::mat>> paramPatterns(
        nParams, std::vector<arma::mat>(npatterns));

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::uword r = paramRow[p];
      const arma::uword c = paramCol[p];

      const arma::vec br = Binv.col(r);
      const arma::vec bc = Binv.col(c);
      arma::mat dVarEta = br * bc.t();

      arma::mat dSigma = L_eta * dVarEta * L_eta.t();

      for (int i = 0; i < npatterns; ++i) {
        const arma::uvec& idx = colidx[i];
        if (idx.is_empty()) {
          paramPatterns[p][i].reset();
          continue;
        }
        if (idx.max() >= dSigma.n_rows) {
          paramPatterns[p][i].reset();
          continue;
        }
        paramPatterns[p][i] = dSigma.submat(idx, idx);
      }
    }

    for (arma::uword p = 0; p < nParams; ++p) {
      for (arma::uword q = p; q < nParams; ++q) {
        double hessContrib = 0.0;

        for (int i = 0; i < npatterns; ++i) {
          if (cache.patternWeights[i] <= DBL_MIN) continue;
          if (cache.sigmaInvPatterns[i].is_empty()) continue;
          const arma::uvec& idx = colidx[i];
          if (idx.is_empty()) continue;

          const arma::mat& Cp = paramPatterns[p][i];
          const arma::mat& Cq = paramPatterns[q][i];
          if (Cp.is_empty() || Cq.is_empty()) continue;

          const arma::mat& A = cache.sigmaInvPatterns[i];
          const arma::vec& diff = cache.diffPatterns[i];
          const arma::mat& S = CovPatterns[j][i];
          if (S.is_empty()) continue;
          const double tg = cache.patternWeights[i];

          arma::mat dA = -A * Cq * A;

          arma::mat dGradSigma =
              0.5 * (dA * S * A + A * S * dA +
                     tg * (dA * diff * diff.t() * A +
                           A * diff * diff.t() * dA) -
                     tg * dA);

          hessContrib += arma::accu(dGradSigma % Cp);
        }

        out.hess(p, q) += hessContrib;
        if (p != q) out.hess(q, p) += hessContrib;
      }
    }
  }

  return out;
}


arma::mat gradPsiComplete(const LMSModel&  M,
                          const arma::mat& V,
                          const std::vector<arma::vec>& TGamma,
                          const std::vector<std::vector<arma::vec>>& MeanPatterns,
                          const std::vector<std::vector<arma::mat>>& CovPatterns,
                          const std::vector<arma::uvec>& colidx,
                          const int npatterns) {
  return psiGradHessComplete(M, V, TGamma, MeanPatterns, CovPatterns,
                             colidx, npatterns).grad;
}


arma::mat gradThetaDeltaComplete(const LMSModel&  M,
                                 const arma::mat& V,
                                 const std::vector<arma::vec>& TGamma,
                                 const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                 const std::vector<std::vector<arma::mat>>& CovPatterns,
                                 const std::vector<arma::uvec>& colidx,
                                 const int npatterns) {
  return thetaDeltaGradHessComplete(M, V, TGamma, MeanPatterns, CovPatterns,
                                    colidx, npatterns).grad;
}


GammaXiDerivatives gammaXiGradHessComplete(const LMSModel&  M,
                                           const arma::mat& V,
                                           const std::vector<arma::vec>& TGamma,
                                           const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                           const std::vector<std::vector<arma::mat>>& CovPatterns,
                                           const std::vector<arma::uvec>& colidx,
                                           const int npatterns) {
  GammaXiDerivatives out;
  out.grad = arma::mat(M.Gx.n_rows, M.Gx.n_cols, arma::fill::zeros);
  out.hess = arma::mat(M.Gx.n_elem, M.Gx.n_elem, arma::fill::zeros);

  const arma::mat Oi = make_Oi(M.k, M.numXis);
  const arma::uword nParams = M.Gx.n_elem;
  std::vector<arma::uword> paramRow(nParams);
  std::vector<arma::uword> paramCol(nParams);
  for (arma::uword c = 0; c < M.Gx.n_cols; ++c) {
    for (arma::uword r = 0; r < M.Gx.n_rows; ++r) {
      const arma::uword idx = r + c * M.Gx.n_rows;
      paramRow[idx] = r;
      paramCol[idx] = c;
    }
  }

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

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const arma::mat& gradSigma = cache.gradSigma;
    if (!cache.hasMeanGrad() && !cache.hasCovGrad()) continue;

    const arma::mat Q   = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat Eta = Binv * Q;
    const arma::uword numXi  = muXi.n_elem;
    const arma::uword numEta = muEta.n_elem;

    std::vector<arma::vec> paramDMu(nParams);
    std::vector<std::vector<arma::mat>> paramSigmaPatterns(
        nParams, std::vector<arma::mat>(npatterns));

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::uword r = paramRow[p];
      const arma::uword c = paramCol[p];

      arma::vec e_r(M.Gx.n_rows, arma::fill::zeros); e_r[r] = 1.0;
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
      out.grad(r, c) += mean_contrib + cov_contrib;

      paramDMu[p] = dMu;
      for (int i = 0; i < npatterns; ++i) {
        const arma::uvec& idx = colidx[i];
        if (idx.is_empty()) {
          paramSigmaPatterns[p][i].reset();
          continue;
        }
        arma::mat sub = dSigma.submat(idx, idx);
        paramSigmaPatterns[p][i] = std::move(sub);
      }
    }

    const bool needMeanHess = cache.hasMeanCurv();
    const bool needCovHess  = cache.hasCovGrad();
    if (!needMeanHess && !needCovHess) continue;

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::vec& dMu_p = paramDMu[p];
      for (arma::uword q = p; q < nParams; ++q) {
        double hessContrib = 0.0;

        if (needMeanHess) {
          const arma::vec& dMu_q = paramDMu[q];
          if (!dMu_p.is_empty() && !dMu_q.is_empty()) {
            const arma::vec curvProd = cache.curvMu * dMu_q;
            hessContrib -= arma::dot(dMu_p, curvProd);
          }
        }

        if (needCovHess) {
          for (int i = 0; i < npatterns; ++i) {
            if (cache.patternWeights[i] <= DBL_MIN) continue;
            if (cache.sigmaInvPatterns[i].is_empty()) continue;
            const arma::mat& Cp = paramSigmaPatterns[p][i];
            if (Cp.is_empty()) continue;
            const arma::uvec& idx = colidx[i];
            if (idx.is_empty()) continue;

            const arma::mat& A = cache.sigmaInvPatterns[i];
            const arma::vec& diff = cache.diffPatterns[i];
            const arma::mat& S = CovPatterns[j][i];
            if (S.is_empty()) continue;
            const double tg = cache.patternWeights[i];

            arma::mat Cq = paramSigmaPatterns[q][i];
            arma::mat dA(A.n_rows, A.n_cols, arma::fill::zeros);
            if (!Cq.is_empty()) dA = -A * Cq * A;

            arma::vec Bq(idx.n_elem, arma::fill::zeros);
            if (!paramDMu[q].is_empty()) {
              if (idx.max() >= paramDMu[q].n_elem) continue;
              Bq = paramDMu[q].elem(idx);
            }
            arma::vec dDiff = -Bq;

            arma::mat dT =
                dA * S * A + A * S * dA +
                tg * (dA * diff * diff.t() * A +
                      A * diff * diff.t() * dA +
                      A * (dDiff * diff.t() + diff * dDiff.t()) * A) -
                tg * dA;

            hessContrib += 0.5 * arma::accu(dT % Cp);
          }
        }

        out.hess(p, q) += hessContrib;
        if (p != q) out.hess(q, p) += hessContrib;
      }
    }
  }

  return out;
}

arma::mat gradGammaXiComplete(const LMSModel&  M,
                              const arma::mat& V,
                              const std::vector<arma::vec>& TGamma,
                              const std::vector<std::vector<arma::vec>>& MeanPatterns,
                              const std::vector<std::vector<arma::mat>>& CovPatterns,
                              const std::vector<arma::uvec>& colidx,
                              const int npatterns) {
  return gammaXiGradHessComplete(M, V, TGamma, MeanPatterns, CovPatterns,
                                 colidx, npatterns).grad;
}


GammaEtaDerivatives gammaEtaGradHessComplete(const LMSModel&  M,
                                             const arma::mat& V,
                                             const std::vector<arma::vec>& TGamma,
                                             const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                             const std::vector<std::vector<arma::mat>>& CovPatterns,
                                             const std::vector<arma::uvec>& colidx,
                                             const int npatterns) {
  GammaEtaDerivatives out;
  out.grad = arma::mat(M.Ge.n_rows, M.Ge.n_cols, arma::fill::zeros);
  out.hess = arma::mat(M.Ge.n_elem, M.Ge.n_elem, arma::fill::zeros);
  const arma::mat Oi = make_Oi(M.k, M.numXis);

  const arma::uword nParams = M.Ge.n_elem;
  std::vector<arma::uword> paramRow(nParams);
  std::vector<arma::uword> paramCol(nParams);
  for (arma::uword c = 0; c < M.Ge.n_cols; ++c) {
    for (arma::uword r = 0; r < M.Ge.n_rows; ++r) {
      const arma::uword idx = r + c * M.Ge.n_rows;
      paramRow[idx] = r;
      paramCol[idx] = c;
    }
  }

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

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const arma::mat& gradSigma = cache.gradSigma;
    if (!cache.hasMeanGrad() && !cache.hasCovGrad()) continue;

    const arma::mat Q   = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat Eta = Binv * Q;
    const arma::uword numXi  = muXi.n_elem;
    const arma::uword numEta = muEta.n_elem;

    std::vector<arma::vec> paramDMu(nParams);
    std::vector<std::vector<arma::mat>> paramSigmaPatterns(
        nParams, std::vector<arma::mat>(npatterns));

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::uword r = paramRow[p];
      const arma::uword c = paramCol[p];
      arma::vec e_r(M.Ge.n_rows, arma::fill::zeros); e_r[r] = 1.0;
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
      out.grad(r, c) += mean_contrib + cov_contrib;

      paramDMu[p] = dMu;
      for (int i = 0; i < npatterns; ++i) {
        const arma::uvec& idx = colidx[i];
        if (idx.is_empty()) {
          paramSigmaPatterns[p][i].reset();
          continue;
        }
        paramSigmaPatterns[p][i] = dSigma.submat(idx, idx);
      }
    }

    const bool needMeanHess = cache.hasMeanCurv();
    const bool needCovHess  = cache.hasCovGrad();
    if (!needMeanHess && !needCovHess) continue;

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::vec& dMu_p = paramDMu[p];
      for (arma::uword q = p; q < nParams; ++q) {
        double hessContrib = 0.0;

        if (needMeanHess) {
          const arma::vec& dMu_q = paramDMu[q];
          if (!dMu_p.is_empty() && !dMu_q.is_empty()) {
            const arma::vec curvProd = cache.curvMu * dMu_q;
            hessContrib -= arma::dot(dMu_p, curvProd);
          }
        }

        if (needCovHess) {
          for (int i = 0; i < npatterns; ++i) {
            if (cache.patternWeights[i] <= DBL_MIN) continue;
            if (cache.sigmaInvPatterns[i].is_empty()) continue;
            const arma::mat& Cp = paramSigmaPatterns[p][i];
            if (Cp.is_empty()) continue;
            const arma::uvec& idx = colidx[i];
            if (idx.is_empty()) continue;

            const arma::mat& A = cache.sigmaInvPatterns[i];
            const arma::vec& diff = cache.diffPatterns[i];
            const arma::mat& S = CovPatterns[j][i];
            if (S.is_empty()) continue;
            const double tg = cache.patternWeights[i];

            const arma::mat& Cq = paramSigmaPatterns[q][i];
            arma::mat dA(A.n_rows, A.n_cols, arma::fill::zeros);
            if (!Cq.is_empty()) dA = -A * Cq * A;

            arma::vec Bq(idx.n_elem, arma::fill::zeros);
            if (!paramDMu[q].is_empty()) {
              if (idx.max() >= paramDMu[q].n_elem) continue;
              Bq = paramDMu[q].elem(idx);
            }
            arma::vec dDiff = -Bq;

            arma::mat dT =
                dA * S * A + A * S * dA +
                tg * (dA * diff * diff.t() * A +
                      A * diff * diff.t() * dA +
                      A * (dDiff * diff.t() + diff * dDiff.t()) * A) -
                tg * dA;

            hessContrib += 0.5 * arma::accu(dT % Cp);
          }
        }

        out.hess(p, q) += hessContrib;
        if (p != q) out.hess(q, p) += hessContrib;
      }
    }
  }

  return out;
}

arma::mat gradGammaEtaComplete(const LMSModel&  M,
                               const arma::mat& V,
                               const std::vector<arma::vec>& TGamma,
                               const std::vector<std::vector<arma::vec>>& MeanPatterns,
                               const std::vector<std::vector<arma::mat>>& CovPatterns,
                               const std::vector<arma::uvec>& colidx,
                               const int npatterns) {
  return gammaEtaGradHessComplete(M, V, TGamma, MeanPatterns, CovPatterns,
                                  colidx, npatterns).grad;
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

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const arma::mat& gradSigma = cache.gradSigma;
    if (!cache.hasMeanGrad() && !cache.hasCovGrad()) continue;

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

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const arma::mat& gradSigma = cache.gradSigma;
    if (!cache.hasMeanGrad() && !cache.hasCovGrad()) continue;

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

  auto read_grad = [&](const arma::mat& G, const arma::uword pos) -> double {
    const arma::uword r = row[pos];
    const arma::uword c = col[pos];
    if (symmetric[pos] && r != c) {
      return G(r, c) + G(c, r);
    }
    return G(r, c);
  };

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
      grad[pos] = read_grad(gradAlpha, pos);
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
      grad[pos] = read_grad(gradPsi, pos);
    }
  }

  // Analytic gradient for thetaDelta block
  if (!idxThetaDelta.is_empty()) {
    const arma::mat gradTheta =
        gradThetaDeltaComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);

    for (std::size_t k = 0; k < idxThetaDelta.n_elem; ++k) {
      const std::size_t pos = idxThetaDelta[k];
      grad[pos] = read_grad(gradTheta, pos);
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
      grad[pos] = read_grad(gradOxx, pos);
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


AlphaDerivatives alphaGradHessComplete(const LMSModel&  M,
                                       const arma::mat& V,
                                       const std::vector<arma::vec>& TGamma,
                                       const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                       const std::vector<std::vector<arma::mat>>& CovPatterns,
                                       const std::vector<arma::uvec>& colidx,
                                       const int npatterns) {
  AlphaDerivatives out;
  out.grad = arma::mat(M.a.n_rows, M.a.n_cols, arma::fill::zeros);
  out.hess = arma::mat(M.a.n_elem, M.a.n_elem, arma::fill::zeros);

  const std::size_t J = V.n_rows;
  const arma::uword numEta = M.Ie.n_rows;
  const arma::uword etaStart = M.numXis;
  const arma::uword etaEnd   = etaStart + numEta - 1;
  const arma::mat L_eta = M.lX.cols(etaStart, etaEnd);

  for (std::size_t j = 0; j < J; ++j) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z    = V.row(j).t();
    const arma::vec zVec = make_zvec(M.k, M.numXis, z);
    const arma::vec muXi = M.beta0 + M.A * zVec;
    const arma::mat kronZ = arma::kron(M.Ie, muXi);
    const arma::mat B    = M.Ie - M.Ge - kronZ.t() * M.Oex;
    const arma::mat Binv = arma::inv(B);
    const arma::vec rhs  = M.a + M.Gx * muXi + kronZ.t() * M.Oxx * muXi;
    const arma::vec muEta = Binv * rhs;

    const arma::vec latent = arma::join_cols(muXi, muEta);
    const arma::vec mu     = M.tX + M.lX * latent;
    const arma::mat Sigma  = M.Sigma(z);

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const arma::mat& curvMu = cache.curvMu;

    if (!cache.hasMeanGrad() && !cache.hasMeanCurv()) continue;

    const arma::vec gradLatent = M.lX.t() * gradMu;
    const arma::vec gradMuEta  = gradLatent.subvec(etaStart, etaEnd);
    const arma::vec gradAlphaNode = Binv.t() * gradMuEta;
    out.grad += arma::reshape(gradAlphaNode, out.grad.n_rows, out.grad.n_cols);

    if (cache.hasMeanCurv()) {
      const arma::mat J = L_eta * Binv;
      arma::mat hessNode = -J.t() * curvMu * J;
      // enforce symmetry to control numeric noise
      hessNode = 0.5 * (hessNode + hessNode.t());
      out.hess += hessNode;
    }
  }

  return out;
}

arma::mat gradAlphaComplete(const LMSModel&  M,
                            const arma::mat& V,
                            const std::vector<arma::vec>& TGamma,
                            const std::vector<std::vector<arma::vec>>& MeanPatterns,
                            const std::vector<std::vector<arma::mat>>& CovPatterns,
                            const std::vector<arma::uvec>& colidx,
                            const int npatterns) {
  return alphaGradHessComplete(M, V, TGamma, MeanPatterns, CovPatterns,
                               colidx, npatterns).grad;
}


Beta0Derivatives beta0GradHessComplete(const LMSModel&  M,
                                       const arma::mat& V,
                                       const std::vector<arma::vec>& TGamma,
                                       const std::vector<std::vector<arma::vec>>& MeanPatterns,
                                       const std::vector<std::vector<arma::mat>>& CovPatterns,
                                       const std::vector<arma::uvec>& colidx,
                                       const int npatterns) {
  Beta0Derivatives out;
  out.grad = arma::mat(M.beta0.n_rows, M.beta0.n_cols, arma::fill::zeros);
  out.hess = arma::mat(M.beta0.n_elem, M.beta0.n_elem, arma::fill::zeros);
  if (M.beta0.is_empty()) return out;

  const arma::mat Oi = make_Oi(M.k, M.numXis);
  const arma::uword nParams = M.beta0.n_elem;
  std::vector<arma::uword> paramRow(nParams);
  std::vector<arma::uword> paramCol(nParams);
  for (arma::uword c = 0; c < M.beta0.n_cols; ++c) {
    for (arma::uword r = 0; r < M.beta0.n_rows; ++r) {
      const arma::uword idx = r + c * M.beta0.n_rows;
      paramRow[idx] = r;
      paramCol[idx] = c;
    }
  }

  const std::size_t J = V.n_rows;

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
    const arma::mat Q = M.Gx * M.A + kronZ.t() * M.Oxx * M.A;
    const arma::mat Eta = Binv * Q;
    const arma::vec latent = arma::join_cols(muXi, muEta);
    const arma::vec mu = M.tX + M.lX * latent;
    const arma::mat Sigma = M.Sigma(z);

    const CompleteCaseDerivatives cache =
        computeCompleteCaseDerivatives(mu, Sigma, TGamma[j],
                                       MeanPatterns[j], CovPatterns[j],
                                       colidx, npatterns);
    const arma::vec& gradMu = cache.gradMu;
    const std::vector<arma::mat>& gradSigmaPatterns = cache.gradSigmaPatterns;
    if (!cache.hasMeanGrad() && !cache.hasCovGrad()) continue;

    const arma::uword lenXi = muXi.n_elem;
    const arma::uword numEta = muEta.n_elem;
    const arma::uword dimXi = lenXi;
    const arma::uword dimEta = numEta;
    const arma::uword dimLatent = dimXi + dimEta;

    std::vector<arma::vec> paramDMu(nParams);
    std::vector<std::vector<arma::mat>> paramSigmaPatterns(
        nParams, std::vector<arma::mat>(npatterns));

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::uword r = paramRow[p];
      const arma::uword c = paramCol[p];
      const arma::uword muPos = r + c * M.beta0.n_rows;
      if (muPos >= lenXi) continue;

      arma::vec d_muXi(lenXi, arma::fill::zeros);
      d_muXi[muPos] = 1.0;

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

      arma::vec dLatent = arma::join_cols(d_muXi, d_muEta);
      arma::vec dMu = M.lX * dLatent;
      arma::mat dSigma = M.lX * dVcov * M.lX.t();

      double cov_contrib = 0.0;
      for (int i = 0; i < npatterns; ++i) {
        if (static_cast<std::size_t>(i) >= gradSigmaPatterns.size()) continue;
        if (gradSigmaPatterns[i].is_empty()) continue;
        const arma::uvec& idx = colidx[i];
        if (idx.is_empty()) continue;
        if (idx.max() >= dSigma.n_rows || idx.max() >= dSigma.n_cols) continue;
        arma::mat dSigmaSub = dSigma.submat(idx, idx);
        cov_contrib += arma::accu(gradSigmaPatterns[i] % dSigmaSub);
        paramSigmaPatterns[p][i] = std::move(dSigmaSub);
      }

      const double mean_contrib = arma::dot(gradMu, dMu);
      out.grad(r, c) += mean_contrib + cov_contrib;
      paramDMu[p] = std::move(dMu);
    }

    const bool needMeanHess = cache.hasMeanCurv();
    const bool needCovHess  = cache.hasCovGrad();
    if (!needMeanHess && !needCovHess) continue;

    for (arma::uword p = 0; p < nParams; ++p) {
      const arma::vec& dMu_p = paramDMu[p];
      for (arma::uword q = p; q < nParams; ++q) {
        double hessContrib = 0.0;

        if (needMeanHess) {
          const arma::vec& dMu_q = paramDMu[q];
          if (!dMu_p.is_empty() && !dMu_q.is_empty()) {
            const arma::vec curvProd = cache.curvMu * dMu_q;
            hessContrib -= arma::dot(dMu_p, curvProd);
          }
        }

        if (needCovHess) {
          for (int i = 0; i < npatterns; ++i) {
            if (cache.patternWeights[i] <= DBL_MIN) continue;
            if (cache.sigmaInvPatterns[i].is_empty()) continue;
            const arma::uvec& idx = colidx[i];
            if (idx.is_empty()) continue;

            const arma::mat& Cp = paramSigmaPatterns[p][i];
            const arma::mat& Cq = paramSigmaPatterns[q][i];
            if (Cp.is_empty() || Cq.is_empty()) continue;

            const arma::mat& A = cache.sigmaInvPatterns[i];
            const arma::vec& diff = cache.diffPatterns[i];
            const arma::mat& S = CovPatterns[j][i];
            if (S.is_empty()) continue;
            const double tg = cache.patternWeights[i];

            arma::mat dA = -A * Cq * A;

            arma::vec Bq(idx.n_elem, arma::fill::zeros);
            if (!paramDMu[q].is_empty()) {
              if (idx.max() >= paramDMu[q].n_elem) continue;
              Bq = paramDMu[q].elem(idx);
            }
            arma::vec dDiff = -Bq;

            arma::mat dGradSigma =
                dA * S * A + A * S * dA +
                tg * (dA * diff * diff.t() * A +
                      A * diff * diff.t() * dA +
                      A * (dDiff * diff.t() + diff * dDiff.t()) * A) -
                tg * dA;
            dGradSigma *= 0.5;

            hessContrib += arma::accu(dGradSigma % Cp);
          }
        }

        out.hess(p, q) += hessContrib;
        if (p != q) out.hess(q, p) += hessContrib;
      }
    }
  }

  return out;
}

arma::mat gradBeta0Complete(const LMSModel&  M,
                            const arma::mat& V,
                            const std::vector<arma::vec>& TGamma,
                            const std::vector<std::vector<arma::vec>>& MeanPatterns,
                            const std::vector<std::vector<arma::mat>>& CovPatterns,
                            const std::vector<arma::uvec>& colidx,
                            const int npatterns) {
  return beta0GradHessComplete(M, V, TGamma, MeanPatterns, CovPatterns,
                               colidx, npatterns).grad;
}


arma::mat gradTauXComplete(const LMSModel&  M,
                           const arma::mat& V,
                           const std::vector<arma::vec>& TGamma,
                           const std::vector<std::vector<arma::vec>>& MeanPatterns,
                           const std::vector<std::vector<arma::mat>>& CovPatterns,
                           const std::vector<arma::uvec>& colidx,
                           const int npatterns) {
  return tauXGradHessComplete(M, V, TGamma, MeanPatterns, CovPatterns,
                              colidx, npatterns).grad;
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


template<class F>
Rcpp::List fdHessCentral(LMSModel&         M,
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
  arma::vec base = getParams(M, block, row, col);
  arma::vec incr =
      arma::max(arma::abs(base), arma::vec(p).fill(minAbsPar)) * relStep;

  arma::vec grad(p, arma::fill::zeros);
  arma::mat Hess(p, p, arma::fill::zeros);

  auto evaluate = [&](const arma::vec& delta) -> double {
    LMSModel Mc = M.thread_clone();
    arma::vec params = base + delta;
    setParams(Mc, block, row, col, symmetric, params);
    return fun(Mc);
  };

  const double f0 = fun(M);

  for (std::size_t i = 0; i < p; ++i) {
    arma::vec delta = arma::zeros<arma::vec>(p);
    delta[i] = incr[i];
    const double f_plus = evaluate(delta);
    delta[i] = -incr[i];
    const double f_minus = evaluate(delta);

    grad[i] = (f_plus - f_minus) / (2.0 * incr[i]);
    Hess(i, i) = (f_plus - 2.0 * f0 + f_minus) / (incr[i] * incr[i]);
  }

  for (std::size_t i = 0; i < p; ++i) {
    for (std::size_t j = i + 1; j < p; ++j) {
      arma::vec delta_pp = arma::zeros<arma::vec>(p);
      delta_pp[i] = incr[i];
      delta_pp[j] = incr[j];
      const double f_pp = evaluate(delta_pp);

      arma::vec delta_pm = delta_pp;
      delta_pm[j] = -incr[j];
      const double f_pm = evaluate(delta_pm);

      arma::vec delta_mp = delta_pp;
      delta_mp[i] = -incr[i];
      const double f_mp = evaluate(delta_mp);

      arma::vec delta_mm = delta_mp;
      delta_mm[j] = -incr[j];
      const double f_mm = evaluate(delta_mm);

      const double hij = (f_pp - f_pm - f_mp + f_mm) /
                         (4.0 * incr[i] * incr[j]);
      Hess(i, j) = hij;
      Hess(j, i) = hij;
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

  Rcpp::List res = fdHessCpp(M, comp_ll, block, row, col, symmetric,
                             relStep, minAbs, ncores);

  arma::vec grad = res["gradient"];
  arma::mat Hess = res["Hessian"];

  const arma::uvec idxAlpha = arma::find(block == 8u);
  if (!idxAlpha.is_empty()) {
    const AlphaDerivatives alphaDerivs =
        alphaGradHessComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);
    const arma::mat& gradAlpha = alphaDerivs.grad;
    const arma::mat& hessAlpha = alphaDerivs.hess;

    const arma::uword alphaRows = M.a.n_rows;
    arma::uvec alphaLocal(idxAlpha.n_elem);

    for (arma::uword k = 0; k < idxAlpha.n_elem; ++k) {
      const arma::uword pos = idxAlpha[k];
      const arma::uword r = row[pos];
      const arma::uword c = col[pos];
      grad[pos] = gradAlpha(r, c);
      alphaLocal[k] = r + c * alphaRows;
    }

    for (arma::uword i = 0; i < idxAlpha.n_elem; ++i) {
      const arma::uword pos_i = idxAlpha[i];
      const arma::uword local_i = alphaLocal[i];
      for (arma::uword j = 0; j < idxAlpha.n_elem; ++j) {
        const arma::uword pos_j = idxAlpha[j];
        const arma::uword local_j = alphaLocal[j];
        Hess(pos_i, pos_j) = hessAlpha(local_i, local_j);
      }
    }
  }

  const arma::uvec idxTauX = arma::find(block == 2u);
  if (!idxTauX.is_empty()) {
    const TauXDerivatives tauDerivs =
        tauXGradHessComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);
    const arma::mat& gradTau = tauDerivs.grad;
    const arma::mat& hessTau = tauDerivs.hess;

    const arma::uword tauRows = M.tX.n_rows;
    arma::uvec tauLocal(idxTauX.n_elem);

    for (arma::uword k = 0; k < idxTauX.n_elem; ++k) {
      const arma::uword pos = idxTauX[k];
      const arma::uword r = row[pos];
      const arma::uword c = col[pos];
      grad[pos] = gradTau(r, c);
      tauLocal[k] = r + c * tauRows;
    }

    for (arma::uword i = 0; i < idxTauX.n_elem; ++i) {
      const arma::uword pos_i = idxTauX[i];
      const arma::uword local_i = tauLocal[i];
      for (arma::uword j = 0; j < idxTauX.n_elem; ++j) {
        const arma::uword pos_j = idxTauX[j];
        const arma::uword local_j = tauLocal[j];
        Hess(pos_i, pos_j) = hessTau(local_i, local_j);
      }
    }
  }

  const arma::uvec idxBeta0 = arma::find(block == 9u);
  if (!idxBeta0.is_empty()) {
    const Beta0Derivatives betaDerivs =
        beta0GradHessComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);
    const arma::mat& gradBeta = betaDerivs.grad;
    const arma::mat& hessBeta = betaDerivs.hess;

    const arma::uword betaRows = M.beta0.n_rows;
    arma::uvec betaLocal(idxBeta0.n_elem);

    for (arma::uword k = 0; k < idxBeta0.n_elem; ++k) {
      const arma::uword pos = idxBeta0[k];
      const arma::uword r = row[pos];
      const arma::uword c = col[pos];
      grad[pos] = gradBeta(r, c);
      betaLocal[k] = r + c * betaRows;
    }

    for (arma::uword i = 0; i < idxBeta0.n_elem; ++i) {
      const arma::uword pos_i = idxBeta0[i];
      const arma::uword local_i = betaLocal[i];
      for (arma::uword j = 0; j < idxBeta0.n_elem; ++j) {
        const arma::uword pos_j = idxBeta0[j];
        const arma::uword local_j = betaLocal[j];
        Hess(pos_i, pos_j) = hessBeta(local_i, local_j);
      }
    }
  }

  const arma::uvec idxLambda = arma::find(block == 0u);
  if (!idxLambda.is_empty()) {
    const LambdaXDerivatives lambdaDerivs =
        lambdaXGradHessComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);
    const arma::mat& gradLambda = lambdaDerivs.grad;
    const arma::mat& hessLambda = lambdaDerivs.hess;

    const arma::uword lambdaRows = M.lX.n_rows;
    arma::uvec lambdaLocal(idxLambda.n_elem);

    for (arma::uword k = 0; k < idxLambda.n_elem; ++k) {
      const arma::uword pos = idxLambda[k];
      const arma::uword r = row[pos];
      const arma::uword c = col[pos];
      grad[pos] = gradLambda(r, c);
      lambdaLocal[k] = r + c * lambdaRows;
    }

    for (arma::uword i = 0; i < idxLambda.n_elem; ++i) {
      const arma::uword pos_i = idxLambda[i];
      const arma::uword local_i = lambdaLocal[i];
      for (arma::uword j = 0; j < idxLambda.n_elem; ++j) {
        const arma::uword pos_j = idxLambda[j];
        const arma::uword local_j = lambdaLocal[j];
        Hess(pos_i, pos_j) = hessLambda(local_i, local_j);
      }
    }
  }

  const arma::uvec idxPsi = arma::find(block == 7u);
  if (!idxPsi.is_empty()) {
    const PsiDerivatives psiDerivs =
        psiGradHessComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);
    const arma::mat& gradPsi = psiDerivs.grad;
    const arma::mat& hessPsi = psiDerivs.hess;

    const arma::uword psiRows = M.Psi.n_rows;
    arma::uvec psiLocal(idxPsi.n_elem);

    for (arma::uword k = 0; k < idxPsi.n_elem; ++k) {
      const arma::uword pos = idxPsi[k];
      const arma::uword r = row[pos];
      const arma::uword c = col[pos];
      grad[pos] = gradPsi(r, c);
      psiLocal[k] = r + c * psiRows;
    }

    for (arma::uword i = 0; i < idxPsi.n_elem; ++i) {
      const arma::uword pos_i = idxPsi[i];
      const arma::uword local_i = psiLocal[i];
      for (arma::uword j = 0; j < idxPsi.n_elem; ++j) {
        const arma::uword pos_j = idxPsi[j];
        const arma::uword local_j = psiLocal[j];
        Hess(pos_i, pos_j) = hessPsi(local_i, local_j);
      }
    }
  }

  const arma::uvec idxTheta = arma::find(block == 4u);
  if (!idxTheta.is_empty()) {
    const ThetaDeltaDerivatives thetaDerivs =
        thetaDeltaGradHessComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);
    const arma::mat& gradTheta = thetaDerivs.grad;
    const arma::mat& hessTheta = thetaDerivs.hess;

    const arma::uword thetaRows = M.d.n_rows;
    arma::uvec thetaLocal(idxTheta.n_elem);

    for (arma::uword k = 0; k < idxTheta.n_elem; ++k) {
      const arma::uword pos = idxTheta[k];
      const arma::uword r = row[pos];
      const arma::uword c = col[pos];
      grad[pos] = gradTheta(r, c);
      thetaLocal[k] = r + c * thetaRows;
    }

    for (arma::uword i = 0; i < idxTheta.n_elem; ++i) {
      const arma::uword pos_i = idxTheta[i];
      const arma::uword local_i = thetaLocal[i];
      for (arma::uword j = 0; j < idxTheta.n_elem; ++j) {
        const arma::uword pos_j = idxTheta[j];
        const arma::uword local_j = thetaLocal[j];
        Hess(pos_i, pos_j) = hessTheta(local_i, local_j);
      }
    }
  }

  const arma::uvec idxGammaXi = arma::find(block == 10u);
  if (!idxGammaXi.is_empty()) {
    const GammaXiDerivatives gammaDerivs =
        gammaXiGradHessComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);
    const arma::mat& gradGx = gammaDerivs.grad;
    const arma::mat& hessGx = gammaDerivs.hess;

    const arma::uword gxRows = M.Gx.n_rows;
    arma::uvec gammaLocal(idxGammaXi.n_elem);

    for (arma::uword k = 0; k < idxGammaXi.n_elem; ++k) {
      const arma::uword pos = idxGammaXi[k];
      const arma::uword r = row[pos];
      const arma::uword c = col[pos];
      grad[pos] = gradGx(r, c);
      gammaLocal[k] = r + c * gxRows;
    }

    for (arma::uword i = 0; i < idxGammaXi.n_elem; ++i) {
      const arma::uword pos_i = idxGammaXi[i];
      const arma::uword local_i = gammaLocal[i];
      for (arma::uword j = 0; j < idxGammaXi.n_elem; ++j) {
        const arma::uword pos_j = idxGammaXi[j];
        const arma::uword local_j = gammaLocal[j];
        Hess(pos_i, pos_j) = hessGx(local_i, local_j);
      }
    }
  }

  const arma::uvec idxGammaEta = arma::find(block == 11u);
  if (!idxGammaEta.is_empty()) {
    const GammaEtaDerivatives gammaDerivs =
        gammaEtaGradHessComplete(M, V, TGamma, Mean, Cov, colidx, npatterns);
    const arma::mat& gradGe = gammaDerivs.grad;
    const arma::mat& hessGe = gammaDerivs.hess;

    const arma::uword geRows = M.Ge.n_rows;
    arma::uvec gammaEtaLocal(idxGammaEta.n_elem);

    for (arma::uword k = 0; k < idxGammaEta.n_elem; ++k) {
      const arma::uword pos = idxGammaEta[k];
      const arma::uword r = row[pos];
      const arma::uword c = col[pos];
      grad[pos] = gradGe(r, c);
      gammaEtaLocal[k] = r + c * geRows;
    }

    for (arma::uword i = 0; i < idxGammaEta.n_elem; ++i) {
      const arma::uword pos_i = idxGammaEta[i];
      const arma::uword local_i = gammaEtaLocal[i];
      for (arma::uword j = 0; j < idxGammaEta.n_elem; ++j) {
        const arma::uword pos_j = idxGammaEta[j];
        const arma::uword local_j = gammaEtaLocal[j];
        Hess(pos_i, pos_j) = hessGe(local_i, local_j);
      }
    }
  }

  res["gradient"] = grad;
  res["Hessian"]  = Hess;
  return res;
}


// [[Rcpp::export]]
Rcpp::List hessCompLogLikLmsCpp_fd(const Rcpp::List& modelR,
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

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(mod, V, TGamma, Mean, Cov,
                                   colidx, n, d, npatterns); // single-threaded
  };

  return fdHessCpp(M, comp_ll, block, row, col, symmetric,
      relStep, minAbs, ncores);
}


// [[Rcpp::export]]
Rcpp::List hessCompLogLikLmsCpp_central(const Rcpp::List& modelR,
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

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(mod, V, TGamma, Mean, Cov,
                                   colidx, n, d, npatterns); // single-threaded
  };

  return fdHessCentral(M, comp_ll, block, row, col, symmetric,
                       relStep, minAbs, ncores);
}
