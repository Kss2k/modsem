#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>

#include "lms.h"
#include "utils.h"
#include "mvnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]


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

  // Forward-pass cache: mu(z)/Sigma(z) together with every intermediate that
  // accumulateMuSigmaAdjointsFromCache() needs to run the reverse pass at the
  // same z. Every hot-loop call site previously called muSigma(z) to get
  // (mu, Sigma), then separately called accumulateMuSigmaAdjoints(M, z, ...),
  // which recomputed this ENTIRE forward pass again from scratch (same z,
  // same model -- an exact, not approximate, duplication of work). Computing
  // it once via muSigmaForward() and threading the cache through removes that
  // duplication; muSigma() below is now a thin wrapper over it so every other
  // caller's behavior/signature is unchanged.
  struct ForwardCache {
    arma::vec zVec, muXi, r, muEta, xieta, mu;
    arma::mat kronZ, B, Binv, C, Eta, varXi, varEta, covXiEta, V, Sigma;
    arma::mat H, R, L; // only populated when hasComposites
  };

  ForwardCache muSigmaForward(const arma::vec& z) const {
    ForwardCache c;
    c.zVec  = makeZvec(k, numXis, z);
    c.muXi  = beta0 + A * c.zVec;
    c.kronZ = arma::kron(Ie, c.muXi);
    c.B     = Ie - Ge - c.kronZ.t() * Oex;

    c.Binv  = arma::inv(c.B);
    c.r     = a + Gx * c.muXi + c.kronZ.t() * Oxx * c.muXi + ZetaProj * c.zVec;
    c.muEta = c.Binv * c.r;
    c.C     = GxA + c.kronZ.t() * Oxx * A + ZetaProj;
    c.Eta   = c.Binv * c.C;
    c.varXi    = AOiAt;
    c.varEta   = c.Eta * Oi * c.Eta.t() + c.Binv * PsiOrth * c.Binv.t();
    c.covXiEta = AOi * c.Eta.t();

    c.V = arma::join_cols(
      arma::join_rows(c.varXi,        c.covXiEta),
      arma::join_rows(c.covXiEta.t(), c.varEta)
    );
    c.xieta = arma::join_cols(c.muXi, c.muEta);

    if (hasComposites) {
      c.H = W.t() * T * W;
      c.R = arma::pinv(c.H);
      c.L = lX + T * W * c.R;
      const arma::mat dc = d + T - c.L * c.H * c.L.t();
      c.mu    = tX + c.L * c.xieta;
      c.Sigma = c.L * c.V * c.L.t() + dc;
    } else {
      c.mu    = tX + lX * c.xieta;
      c.Sigma = lX * c.V * lX.t() + d;
    }

    return c;
  }

  // Combined mu+Sigma: avoids recomputing zVec/kronZ/Binv/lXc for the same z.
  // Use this in any hot path that needs both but not the reverse pass.
  std::pair<arma::vec, arma::mat> muSigma(const arma::vec& z) const {
    const ForwardCache c = muSigmaForward(z);
    return std::make_pair(c.mu, c.Sigma);
  }

  // mu(z), Sigma(z) together with their first derivatives w.r.t. z (k slices
  // each). mu(z) = tX + L*xieta(z), Sigma(z) = L*V(z)*L' + dConst, where L
  // (lX or, with composites, the composite loading lXc) and dConst are
  // z-independent, so only xieta(z) and V(z) need differentiating. B(z) is
  // affine in z (dB/dz_l is a constant matrix), so d(Binv)/dz_l follows the
  // standard matrix-inverse identity -Binv*dB_l*Binv. Mirrors muSigma() term
  // for term, just forward-differentiated w.r.t. each z_l.
  struct MuSigmaJacZ {
    arma::vec mu;
    arma::mat Sigma;
    std::vector<arma::vec> dmu;    // length k
    std::vector<arma::mat> dSigma; // length k
  };

  MuSigmaJacZ muSigmaJacZ(const arma::vec& z) const {
    const arma::vec zVec  = makeZvec(k, numXis, z);
    const arma::vec muXi  = beta0 + A * zVec;
    const arma::mat kronZ = arma::kron(Ie, muXi);
    const arma::mat B     = Ie - Ge - kronZ.t() * Oex;

    const arma::mat Binv  = arma::inv(B);
    const arma::vec r     = a + Gx * muXi + kronZ.t() * Oxx * muXi + ZetaProj * zVec;
    const arma::vec muEta = Binv * r;
    const arma::mat C     = GxA + kronZ.t() * Oxx * A + ZetaProj;
    const arma::mat Eta   = Binv * C;
    const arma::mat varXi    = AOiAt;
    const arma::mat varEta   = Eta * Oi * Eta.t() + Binv * PsiOrth * Binv.t();
    const arma::mat covXiEta = AOi * Eta.t();

    const arma::mat V = arma::join_cols(
      arma::join_rows(varXi,        covXiEta),
      arma::join_rows(covXiEta.t(), varEta)
    );
    const arma::vec xieta = arma::join_cols(muXi, muEta);

    arma::mat L, dConst;
    if (hasComposites) {
      const arma::mat H = W.t() * T * W;
      const arma::mat R = arma::pinv(H);
      L = lX + T * W * R;
      dConst = d + T - L * W.t() * T * W * L.t();
    } else {
      L = lX;
      dConst = d;
    }

    MuSigmaJacZ out;
    out.mu    = tX + L * xieta;
    out.Sigma = L * V * L.t() + dConst;
    out.dmu.resize(k);
    out.dSigma.resize(k);

    for (unsigned l = 0; l < k; ++l) {
      const arma::vec dMuXi_l = A.col(l);              // d(muXi)/dz_l
      const arma::mat K_l     = arma::kron(Ie, dMuXi_l); // d(kronZ)/dz_l
      const arma::mat B_l     = -K_l.t() * Oex;          // d(B)/dz_l, constant in z
      const arma::mat dBinv_l = -Binv * B_l * Binv;

      const arma::vec dr_l = Gx * dMuXi_l + K_l.t() * Oxx * muXi
                            + kronZ.t() * Oxx * dMuXi_l + ZetaProj.col(l);
      const arma::vec dMuEta_l = dBinv_l * r + Binv * dr_l;

      const arma::mat dC_l   = K_l.t() * Oxx * A;
      const arma::mat dEta_l = dBinv_l * C + Binv * dC_l;

      const arma::mat dCovXiEta_l = AOi * dEta_l.t();
      const arma::mat dVarEta_l   = dEta_l * Oi * Eta.t() + Eta * Oi * dEta_l.t()
                                   + dBinv_l * PsiOrth * Binv.t() + Binv * PsiOrth * dBinv_l.t();

      const arma::mat dVarXi_l(varXi.n_rows, varXi.n_cols, arma::fill::zeros);
      const arma::mat dV_l = arma::join_cols(
        arma::join_rows(dVarXi_l,        dCovXiEta_l),
        arma::join_rows(dCovXiEta_l.t(), dVarEta_l)
      );

      const arma::vec dxieta_l = arma::join_cols(dMuXi_l, dMuEta_l);

      out.dmu[l]    = L * dxieta_l;
      out.dSigma[l] = L * dV_l * L.t();
    }

    return out;
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


// Evaluate exp(MVN density) for every row of X under N(mu[idx], Sig[idx,idx]).
// No OMP thread-count manipulation — safe to call from inside a parallel region.
static arma::vec mvnDensSt(const arma::mat&  X,
                            const arma::uvec& idx,
                            const arma::vec&  mu,
                            const arma::mat&  Sig) {
  const arma::vec muS  = mu.elem(idx);
  const arma::mat sigS = Sig.submat(idx, idx);

  arma::mat L;
  if (!sigS.is_finite() || !arma::chol(L, sigS, "lower")) { // avoid warning in arma::chol
    arma::vec out(X.n_rows);
    out.fill(arma::datum::nan);
    return out;
  }

  const arma::uword d  = idx.n_elem;
  const double logDet  = 2.0 * arma::sum(arma::log(L.diag()));
  const double c       = -0.5 * (static_cast<double>(d) * std::log(2.0 * M_PI) + logDet);

  // Solve L * Z^T = (X - muS)^T via forward substitution — avoids row-by-row loop.
  const arma::mat cent = X.each_row() - muS.t();
  const arma::mat Zt   = arma::solve(arma::trimatl(L), cent.t(), arma::solve_opts::fast);

  return arma::exp(c - 0.5 * arma::sum(arma::square(Zt), 0).t());
}


// [[Rcpp::export]]
Rcpp::List muSigmaLmsCpp(Rcpp::List model, arma::vec z) {
  const std::pair<arma::vec, arma::mat> ms = LMSModel(model).muSigma(z);
  return Rcpp::List::create(
    Rcpp::Named("mu")    = ms.first,
    Rcpp::Named("sigma") = ms.second
  );
}


// Gradient of log f(y | z) w.r.t. z (dimension k), for ONE observation, using
// the analytic Jacobians from muSigmaJacZ(). `y`/`idx` are the observed sub-
// vector and its column indices (missing-data pattern). Does NOT include the
// -0.5*z'z prior term (added separately, trivially, wherever this is used).
// Returns a k-vector of NaN if Sigma restricted to idx is not PD.
inline arma::vec gradLogDensityZ(const LMSModel& M,
                                 const arma::vec& z,
                                 const arma::vec& y,
                                 const arma::uvec& idx) {
  const LMSModel::MuSigmaJacZ jac = M.muSigmaJacZ(z);
  const arma::vec muS  = jac.mu.elem(idx);
  const arma::mat sigS = jac.Sigma.submat(idx, idx);

  arma::vec grad(M.k);

  arma::mat L;
  if (!sigS.is_finite() || !arma::chol(L, sigS, "lower")) {
    grad.fill(arma::datum::nan);
    return grad;
  }

  const arma::mat Sinv = arma::solve(
    arma::trimatu(L.t()),
    arma::solve(arma::trimatl(L),
                arma::eye<arma::mat>(sigS.n_rows, sigS.n_cols),
                arma::solve_opts::fast),
    arma::solve_opts::fast
  );
  const arma::vec diff     = y - muS;
  const arma::vec Sinv_diff = Sinv * diff;

  for (unsigned l = 0; l < M.k; ++l) {
    const arma::vec dmuS  = jac.dmu[l].elem(idx);
    const arma::mat dSigS = jac.dSigma[l].submat(idx, idx);

    const double trTerm   = arma::trace(Sinv * dSigS);
    const double quadTerm = arma::as_scalar(Sinv_diff.t() * dSigS * Sinv_diff);

    grad[l] = -0.5 * trTerm + arma::as_scalar(diff.t() * Sinv * dmuS) + 0.5 * quadTerm;
  }

  return grad;
}


// [[Rcpp::export]]
arma::vec gradLogDensityZLmsCpp(const Rcpp::List& modelR, const arma::vec& z,
                                const arma::vec& y, const arma::uvec& idx) {
  const LMSModel M(modelR);
  return gradLogDensityZ(M, z, y, idx - 1); // idx passed 1-indexed from R
}


// ======================================================
// Per-row AGHQ: Newton mode-finding + adaptive node/weight construction
// ======================================================

struct AghqRowResult {
  arma::vec z;   // mode (or fallback 0)
  arma::mat L;   // lower Cholesky of curvature Hneg, Hneg = L*L' (or fallback I)
  bool      ok;  // whether Newton actually converged (false => fallback used)
};


// Newton mode-search for g(z) = log f(y|z) - 0.5*z'z, using the analytic
// gradient (gradLogDensityZ) and a Hessian obtained by central-differencing
// that analytic gradient (2k evaluations per iteration; k is small so this
// is cheap). Falls back to z*=0, curvature=I (i.e. the same result as
// "fixed-n") on any failure: non-finite density, a Hessian that stays
// indefinite even after ridging, or hitting `iterMax` without reaching
// `gradTol`.
inline AghqRowResult findModeZ(const LMSModel&   M,
                               const arma::vec&  y,
                               const arma::uvec& idx,
                               const arma::vec&  z0,
                               const unsigned    iterMax = 25,
                               const double      gradTol = 1e-8,
                               const double      fdH     = 1e-4) {
  const unsigned k = M.k;

  AghqRowResult fallback;
  fallback.z  = arma::zeros<arma::vec>(k);
  fallback.L  = arma::eye<arma::mat>(k, k);
  fallback.ok = false;

  if (k == 0) return fallback;

  auto gradG = [&](const arma::vec& zz) -> arma::vec {
    arma::vec g = gradLogDensityZ(M, zz, y, idx);
    if (!g.is_finite()) return g;
    return g - zz; // add d(-0.5 z'z)/dz = -z
  };

  auto hessNeg = [&](const arma::vec& zz, const arma::vec& gCurrent,
                     arma::mat& Hneg) -> bool {
    Hneg.zeros(k, k);
    for (unsigned j = 0; j < k; ++j) {
      arma::vec zp = zz, zm = zz;
      zp[j] += fdH; zm[j] -= fdH;
      const arma::vec gp = gradG(zp);
      const arma::vec gm = gradG(zm);
      if (!gp.is_finite() || !gm.is_finite()) return false;
      Hneg.col(j) = -(gp - gm) / (2.0 * fdH);
    }
    Hneg = 0.5 * (Hneg + Hneg.t());
    return true;
  };

  auto cholWithRidge = [&](arma::mat& Hneg, arma::mat& L) -> bool {
    double ridge = 0.0;
    for (int attempt = 0; attempt < 6; ++attempt) {
      arma::mat Htry = Hneg + ridge * arma::eye<arma::mat>(k, k);
      if (arma::chol(L, Htry, "lower")) { Hneg = Htry; return true; }
      ridge = (ridge <= 0.0) ? 1e-6 : ridge * 10.0;
    }
    return false;
  };

  arma::vec z = (z0.n_elem == k && z0.is_finite()) ? z0 : arma::zeros<arma::vec>(k);

  for (unsigned it = 0; it < iterMax; ++it) {
    const arma::vec g = gradG(z);
    if (!g.is_finite()) return fallback;

    arma::mat Hneg;
    if (!hessNeg(z, g, Hneg)) return fallback;

    arma::mat L;
    if (!cholWithRidge(Hneg, L)) return fallback;

    if (arma::abs(g).max() < gradTol) {
      AghqRowResult out; out.z = z; out.L = L; out.ok = true;
      return out;
    }

    // Newton step: Hneg * delta = g
    arma::vec delta = arma::solve(
      arma::trimatu(L.t()),
      arma::solve(arma::trimatl(L), g, arma::solve_opts::fast),
      arma::solve_opts::fast
    );

    double step = 1.0;
    arma::vec zNew = z + delta;
    arma::vec gNew = gradG(zNew);
    int halvings = 0;
    const double gNorm = arma::abs(g).max();
    while ((!gNew.is_finite() || arma::abs(gNew).max() > 10.0 * gNorm + 1.0) &&
           halvings < 8) {
      step *= 0.5;
      zNew = z + step * delta;
      gNew = gradG(zNew);
      ++halvings;
    }
    if (!gNew.is_finite()) return fallback;

    z = zNew;
  }

  return fallback; // iteration cap reached without meeting gradTol
}


struct AghqGrid {
  arma::mat Z; // J x k adapted nodes
  arma::vec W; // J adapted weights
};


// Builds one row's adapted nodes/weights from mode `zStar` and curvature
// Cholesky `L` (Hneg = L*L'), given the shared base grid (baseN, baseW) --
// the SAME nodes/weights `quadrature()` produces on the R side (i.e. already
// sqrt(2)-scaled, pi^(-k/2)-weighted). With zStar=0, L=I this reduces
// exactly to (baseN, baseW): see derivation notes in the R-side aghq wiring.
inline AghqGrid buildAghqGrid(const arma::vec& zStar, const arma::mat& L,
                              const arma::mat& baseN, const arma::vec& baseW) {
  const unsigned k = zStar.n_elem;
  const unsigned J = baseN.n_rows;

  const arma::mat Linv = arma::inv(arma::trimatl(L)); // L^{-1}
  const arma::mat Z = arma::repmat(zStar.t(), J, 1) + baseN * Linv;

  const double logDetTerm = k * 0.5 * std::log(2.0 * M_PI) - arma::sum(arma::log(L.diag()));

  arma::vec W(J);
  for (unsigned j = 0; j < J; ++j) {
    const double nsq   = arma::dot(baseN.row(j), baseN.row(j));
    const double zsq   = arma::dot(Z.row(j), Z.row(j));
    const double logphi = -0.5 * (k * std::log(2.0 * M_PI) + zsq);
    W[j] = std::exp(logDetTerm + std::log(baseW[j]) + 0.5 * nsq + logphi);
  }

  AghqGrid out; out.Z = Z; out.W = W;
  return out;
}


// For every row (grouped by missing-data pattern), finds the AGHQ mode +
// curvature (or, if `adapt=false`, uses the fixed-n fallback z*=0, I for
// every row) and builds that row's adapted node/weight grid. `Z0` holds
// warm-start modes in the same row order as the concatenated patterns
// (i.e. the order data$ids/data$data.split are iterated in elsewhere).
// [[Rcpp::export]]
Rcpp::List aghqRowGridsLmsCpp(const Rcpp::List& modelR,
                              const Rcpp::List& dataR,
                              const Rcpp::List& colidxR,
                              const arma::uvec& n,
                              const arma::mat&  baseN,
                              const arma::vec&  baseW,
                              const arma::mat&  Z0,
                              const bool         adapt,
                              const int          npatterns = 1,
                              const unsigned      iterMax  = 25,
                              const double        gradTol  = 1e-8,
                              const double        fdH      = 1e-4,
                              const int          ncores    = 1) {
  const LMSModel M(modelR);
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);
  const unsigned k  = M.k;
  const unsigned J  = baseN.n_rows;
  const arma::uword N = arma::sum(n);

  arma::cube Z(N, J, k, arma::fill::zeros);
  arma::mat  W(N, J, arma::fill::zeros);
  arma::mat  Zstar(N, k, arma::fill::zeros);
  arma::uvec ok(N, arma::fill::zeros);

  ThreadSetter ts(ncores);

  arma::uword offset = 0;
  for (int pat = 0; pat < npatterns; ++pat) {
    const arma::uvec& idx = colidx[pat];
    const arma::mat&  Dat = data[pat];
    const arma::uword nPat = n[pat];

#pragma omp parallel for schedule(static) default(none) if(ncores > 1) \
    shared(M, Dat, idx, baseN, baseW, Z0, adapt, offset, nPat, k, J, \
           Z, W, Zstar, ok, iterMax, gradTol, fdH)
    for (arma::uword i = 0; i < nPat; ++i) {
      const arma::uword row = offset + i;
      const arma::vec y  = Dat.row(i).t();
      const arma::vec z0 = Z0.row(row).t();

      AghqRowResult res;
      if (adapt) {
        res = findModeZ(M, y, idx, z0, iterMax, gradTol, fdH);
      } else {
        res.z  = arma::zeros<arma::vec>(k);
        res.L  = arma::eye<arma::mat>(k, k);
        res.ok = true;
      }

      const AghqGrid grid = buildAghqGrid(res.z, res.L, baseN, baseW);

      for (arma::uword j = 0; j < J; ++j)
        for (arma::uword l = 0; l < k; ++l)
          Z(row, j, l) = grid.Z(j, l);

      W.row(row) = grid.W.t();
      Zstar.row(row) = res.z.t();
      ok[row] = res.ok ? 1 : 0;
    }
    offset += nPat;
  }

  return Rcpp::List::create(
    Rcpp::Named("Z")     = Z,
    Rcpp::Named("W")     = W,
    Rcpp::Named("Zstar") = Zstar,
    Rcpp::Named("ok")    = ok
  );
}


inline double& lmsParam(LMSModel& M, std::size_t blk,
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
    case 14: return  M.W    (r,c);
    case 15: return  M.T    (r,c);
    case 17: return  M.covZetaXi (r,c);
    default: Rcpp::stop("unknown block id");
  }
}

struct LMSAdjoints {
  arma::mat A, Oxx, Oex, lX, tX, Gx, Ge, a, beta0, Psi, covZetaXi, d, W, T;

  explicit LMSAdjoints(const LMSModel& M) :
    A(M.A.n_rows, M.A.n_cols, arma::fill::zeros),
    Oxx(M.Oxx.n_rows, M.Oxx.n_cols, arma::fill::zeros),
    Oex(M.Oex.n_rows, M.Oex.n_cols, arma::fill::zeros),
    lX(M.lX.n_rows, M.lX.n_cols, arma::fill::zeros),
    tX(M.tX.n_rows, M.tX.n_cols, arma::fill::zeros),
    Gx(M.Gx.n_rows, M.Gx.n_cols, arma::fill::zeros),
    Ge(M.Ge.n_rows, M.Ge.n_cols, arma::fill::zeros),
    a(M.a.n_rows, M.a.n_cols, arma::fill::zeros),
    beta0(M.beta0.n_rows, M.beta0.n_cols, arma::fill::zeros),
    Psi(M.Psi.n_rows, M.Psi.n_cols, arma::fill::zeros),
    covZetaXi(M.covZetaXi.n_rows, M.covZetaXi.n_cols, arma::fill::zeros),
    d(M.d.n_rows, M.d.n_cols, arma::fill::zeros),
    W(M.W.n_rows, M.W.n_cols, arma::fill::zeros),
    T(M.T.n_rows, M.T.n_cols, arma::fill::zeros) {}

  void add(const LMSAdjoints& other) {
    A         += other.A;
    Oxx       += other.Oxx;
    Oex       += other.Oex;
    lX        += other.lX;
    tX        += other.tX;
    Gx        += other.Gx;
    Ge        += other.Ge;
    a         += other.a;
    beta0     += other.beta0;
    Psi       += other.Psi;
    covZetaXi += other.covZetaXi;
    d         += other.d;
    W         += other.W;
    T         += other.T;
  }
};


inline void addKronIeUAdjoint(arma::vec& uBar,
                              const arma::mat& KBar,
                              const unsigned numXis,
                              const unsigned numEtas) {
  for (unsigned e = 0; e < numEtas; ++e) {
    uBar += KBar.submat(e * numXis, e, (e + 1) * numXis - 1, e);
  }
}


inline const arma::mat& lmsAdjointBlock(const LMSAdjoints& adj,
                                        std::size_t blk) {
  switch (blk) {
    case 0 : return adj.lX;
    case 2 : return adj.tX;
    case 4 : return adj.d;
    case 6 : return adj.A;
    case 7 : return adj.Psi;
    case 8 : return adj.a;
    case 9 : return adj.beta0;
    case 10: return adj.Gx;
    case 11: return adj.Ge;
    case 12: return adj.Oxx;
    case 13: return adj.Oex;
    case 14: return adj.W;
    case 15: return adj.T;
    case 17: return adj.covZetaXi;
    default: {
      static const arma::mat empty;
      return empty;
    }
  }
}


inline void addLatentCovAdjoints(const LMSModel& M,
                                 const arma::mat& ZetaProjBar,
                                 const arma::mat& psiIndependentBar,
                                 LMSAdjoints& adj) {
  arma::mat EBar = ZetaProjBar;

  if (psiIndependentBar.n_elem) {
    adj.Psi += psiIndependentBar;
    EBar -= (psiIndependentBar + psiIndependentBar.t()) * M.ZetaProj;
  }

  if (!EBar.n_elem) return;

  const arma::mat invA = arma::inv(arma::trimatl(M.A));
  adj.covZetaXi += EBar * invA;
  adj.A    -= M.ZetaProj.t() * EBar * invA;
}


inline void accumulateMuSigmaAdjointsFromCache(const LMSModel& M,
                                               const LMSModel::ForwardCache& c,
                                               const arma::vec& muBar,
                                               const arma::mat& SigmaBar,
                                               LMSAdjoints& adj) {
  const arma::uword numEta = M.Ie.n_rows;
  const arma::vec& zVec = c.zVec;
  const arma::vec& u    = c.muXi;
  const arma::mat& K    = c.kronZ;
  const arma::mat& Binv = c.Binv;
  const arma::vec& r    = c.r;
  const arma::mat& C    = c.C;
  const arma::mat& Eta  = c.Eta;
  const arma::mat& V    = c.V;
  const arma::vec& xieta = c.xieta;

  if (adj.tX.n_elem == muBar.n_elem) adj.tX += muBar;
  arma::vec xietaBar;
  arma::mat VBar;

  if (M.hasComposites) {
    // mu = L*xieta and Sigma = L*(V-H)*L' + d + T, where
    // H = W'*T*W and L = lX + T*W*pinv(H).
    const arma::mat& H = c.H;
    const arma::mat& R = c.R;
    const arma::mat& L = c.L;
    const arma::mat VC = V - H;

    arma::mat LBar = muBar * xieta.t() +
      (SigmaBar + SigmaBar.t()) * L * VC;
    arma::mat HBar = -L.t() * SigmaBar * L;
    const arma::mat RBar = M.W.t() * M.T.t() * LBar;

    adj.lX += LBar;
    adj.d  += SigmaBar;
    adj.T  += SigmaBar + LBar * R.t() * M.W.t();
    adj.W  += M.T.t() * LBar * R.t();

    // Reverse the Moore-Penrose inverse. This is the constant-rank
    // differential and reduces to -R' RBar R' when H is nonsingular.
    const arma::mat IH = arma::eye<arma::mat>(H.n_rows, H.n_cols);
    const arma::mat leftNull  = IH - R * H;
    const arma::mat rightNull = IH - H * R;
    HBar += -R.t() * RBar * R.t()
      + (R * R.t()) * RBar * rightNull.t()
      + leftNull.t() * RBar * (R.t() * R);

    // H = W'*T*W.
    adj.T += M.W * HBar * M.W.t();
    adj.W += M.T * M.W * HBar.t() + M.T.t() * M.W * HBar;

    xietaBar = L.t() * muBar;
    VBar = L.t() * SigmaBar * L;
  } else {
    adj.lX += muBar * xieta.t();
    adj.lX += (SigmaBar + SigmaBar.t()) * M.lX * V;
    adj.d  += SigmaBar;
    xietaBar = M.lX.t() * muBar;
    VBar = M.lX.t() * SigmaBar * M.lX;
  }

  arma::vec uBar = xietaBar.subvec(0, M.numXis - 1);
  arma::vec vBar = xietaBar.subvec(M.numXis, M.numXis + numEta - 1);

  arma::mat varXiBar = VBar.submat(0, 0, M.numXis - 1, M.numXis - 1);
  arma::mat covBar = VBar.submat(0, M.numXis, M.numXis - 1, M.numXis + numEta - 1);
  arma::mat varEtaBar = VBar.submat(M.numXis, M.numXis,
                                    M.numXis + numEta - 1,
                                    M.numXis + numEta - 1);
  covBar += VBar.submat(M.numXis, 0,
                        M.numXis + numEta - 1,
                        M.numXis - 1).t();

  adj.A += (varXiBar + varXiBar.t()) * M.AOi;

  arma::mat EtaBar = covBar.t() * M.AOi;
  adj.A += covBar * Eta * M.Oi;

  EtaBar += (varEtaBar + varEtaBar.t()) * Eta * M.Oi;
  arma::mat BinvBar = varEtaBar * Binv * M.PsiOrth.t() +
                      varEtaBar.t() * Binv * M.PsiOrth;
  arma::mat PsiOrthBar = Binv.t() * varEtaBar * Binv;

  BinvBar += EtaBar * C.t();
  arma::mat CBar = Binv.t() * EtaBar;

  BinvBar += vBar * r.t();
  arma::vec rBar = Binv.t() * vBar;

  arma::mat BBar = -Binv.t() * BinvBar * Binv.t();
  adj.Ge -= BBar;
  arma::mat KBar = -M.Oex * BBar.t();
  adj.Oex += -K * BBar;

  adj.a += rBar;
  adj.Gx += rBar * u.t();
  uBar += M.Gx.t() * rBar;
  arma::mat ZetaProjBar = rBar * zVec.t();

  arma::vec q = M.Oxx * u;
  KBar += q * rBar.t();
  arma::vec qBar = K * rBar;
  adj.Oxx += qBar * u.t();
  uBar += M.Oxx.t() * qBar;

  adj.Gx += CBar * M.A.t();
  adj.A += M.Gx.t() * CBar;
  ZetaProjBar += CBar;

  arma::mat MBar = K * CBar;
  KBar += (M.Oxx * M.A) * CBar.t();
  adj.Oxx += MBar * M.A.t();
  adj.A += M.Oxx.t() * MBar;

  addKronIeUAdjoint(uBar, KBar, M.numXis, numEta);
  addLatentCovAdjoints(M, ZetaProjBar, PsiOrthBar, adj);

  adj.beta0 += uBar;
  adj.A += uBar * zVec.t();
}


inline bool completeLogLikScore(const arma::vec& mu,
                                const arma::mat& sigma,
                                const arma::vec& nu,
                                const arma::mat& S,
                                const double tgamma,
                                arma::vec& muBar,
                                arma::mat& SigmaBar) {
  if (!sigma.is_finite()) return false;

  arma::mat L;
  if (!arma::chol(L, sigma, "lower")) return false;

  const arma::vec diff = nu - mu;
  const arma::vec invDiff = arma::solve(
    arma::trimatu(L.t()),
    arma::solve(arma::trimatl(L), diff, arma::solve_opts::fast),
    arma::solve_opts::fast
  );
  const arma::mat Sinv = arma::solve(
    arma::trimatu(L.t()),
    arma::solve(arma::trimatl(L),
                arma::eye<arma::mat>(sigma.n_rows, sigma.n_cols),
                arma::solve_opts::fast),
    arma::solve_opts::fast
  );

  muBar    = tgamma * invDiff;
  SigmaBar = 0.5 * (Sinv * (S + tgamma * diff * diff.t()) * Sinv -
                    tgamma * Sinv);

  return true;
}


// ======================================================
// Per-row AGHQ E-step (posterior) and M-step (complete log-lik/gradient)
// ======================================================
//
// Every row has its own adapted node/weight grid (Z: N x J x k, from
// aghqRowGridsLmsCpp), so unlike the shared-grid path there is no pooling
// across rows at a node: each (row, node) pair is a "group" of exactly one
// observation. The weighted-Gaussian complete-data log-lik/score machinery
// (totalDmvnWeighted / completeLogLikScore) already handles this correctly
// when passed nu=y_i, S=0, tg=gamma_ij -- the general weighted formula
// collapses to the plain single-point case, so it's reused unchanged; no
// sufficient-statistic aggregation (estepSuffStatLmsCpp's analog) is needed.


// E-step: posterior gamma_ij (N x J) from per-row adapted nodes/weights, plus
// the observed (marginal) log-likelihood. Sampling weights, if present, are
// folded into gamma (as in estepLmsGroup) but must not alter within-row
// normalization.
// [[Rcpp::export]]
Rcpp::List aghqEstepLmsCpp(const Rcpp::List& modelR,
                           const arma::cube& Z,
                           const arma::mat&  W,
                           const Rcpp::List& dataR,
                           const Rcpp::List& colidxR,
                           const arma::vec&  samplingWeights,
                           const arma::uvec& n,
                           const int npatterns = 1,
                           const int ncores    = 1) {
  const LMSModel M(modelR);
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);
  const arma::uword J = Z.n_cols;
  const arma::uword k = Z.n_slices;
  const arma::uword N = arma::sum(n);
  const bool hasSW = samplingWeights.n_elem == N;

  arma::mat Gamma(N, J, arma::fill::zeros);
  arma::vec density(N, arma::fill::zeros);

  ThreadSetter ts(ncores);

  arma::uword offset = 0;
  for (int pat = 0; pat < npatterns; ++pat) {
    const arma::uvec& idx  = colidx[pat];
    const arma::mat&  Dat  = data[pat];
    const arma::uword nPat = n[pat];

#pragma omp parallel for schedule(static) default(none) if(ncores > 1) \
    shared(M, Dat, idx, Z, W, offset, nPat, J, k, Gamma, density)
    for (arma::uword i = 0; i < nPat; ++i) {
      const arma::uword row = offset + i;

      arma::vec dens(J);
      for (arma::uword j = 0; j < J; ++j) {
        arma::vec zj(k);
        for (arma::uword l = 0; l < k; ++l) zj[l] = Z(row, j, l);

        const auto ms = M.muSigma(zj);
        // mvnDensSt evaluates the (non-log) MVN density for every row of its
        // first argument; called here with a single row.
        const arma::vec d1 = mvnDensSt(Dat.rows(i, i), idx, ms.first, ms.second);
        const double dj = d1[0];
        dens[j] = std::isfinite(dj) ? W(row, j) * dj : 0.0;
      }

      const double rowDensity = arma::sum(dens);
      density[row] = rowDensity;
      if (rowDensity > DBL_MIN) Gamma.row(row) = (dens / rowDensity).t();
    }
    offset += nPat;
  }

  const arma::vec densitySafe = arma::clamp(density, DBL_MIN, arma::datum::inf);
  const double obsLL = hasSW
    ? arma::sum(samplingWeights % arma::log(densitySafe))
    : arma::sum(arma::log(densitySafe));

  if (hasSW) Gamma.each_col() %= samplingWeights;

  return Rcpp::List::create(
    Rcpp::Named("Gamma")   = Gamma,
    Rcpp::Named("density") = density,
    Rcpp::Named("obsLL")   = obsLL
  );
}


// M-step objective: complete-data log-lik, sum_i sum_j gamma_ij * log f(y_i|z_ij).
// Internal (LMSModel-taking) version, reused by the Rcpp-exported wrapper
// below and by the Hessian-via-FD wrapper further down, to avoid re-parsing
// modelR/dataR/colidxR from R on every FD evaluation.
inline double completeLogLikFromModelAghq(const LMSModel&  M,
                                          const arma::cube& Z,
                                          const arma::mat&  Gamma,
                                          const std::vector<arma::mat>&  data,
                                          const std::vector<arma::uvec>& colidx,
                                          const arma::uvec& n,
                                          const int npatterns = 1,
                                          const int ncores    = 1) {
  const arma::uword J = Z.n_cols;
  const arma::uword k = Z.n_slices;

  ThreadSetter ts(ncores);

  double ll = 0.0;
  arma::uword offset = 0;

  for (int pat = 0; pat < npatterns; ++pat) {
    const arma::uvec& idx  = colidx[pat];
    const arma::mat&  Dat  = data[pat];
    const arma::uword nPat = n[pat];
    const int dObs = (int) idx.n_elem;
    const arma::mat Szero  = arma::zeros<arma::mat>(dObs, dObs);

    double llPat = 0.0;
#pragma omp parallel for reduction(+:llPat) schedule(static) default(none) \
    shared(M, Dat, idx, Z, Gamma, offset, nPat, J, k, Szero, dObs) if(ncores > 1)
    for (arma::uword i = 0; i < nPat; ++i) {
      const arma::uword row = offset + i;
      const arma::vec y = Dat.row(i).t();
      double llRow = 0.0;

      for (arma::uword j = 0; j < J; ++j) {
        const double tg = Gamma(row, j);
        if (tg <= DBL_MIN) continue;

        arma::vec zj(k);
        for (arma::uword l = 0; l < k; ++l) zj[l] = Z(row, j, l);

        const auto ms = M.muSigma(zj);
        const arma::vec muS  = ms.first.elem(idx);
        const arma::mat sigS = ms.second.submat(idx, idx);

        llRow += totalDmvnWeighted(muS, sigS, y, Szero, tg, dObs);
      }
      llPat += llRow;
    }

    ll += llPat;
    offset += nPat;
  }

  return ll;
}


// [[Rcpp::export]]
double completeLogLikLmsAghqCpp(const Rcpp::List& modelR,
                                const arma::cube& Z,
                                const arma::mat&  Gamma,
                                const Rcpp::List& dataR,
                                const Rcpp::List& colidxR,
                                const arma::uvec& n,
                                const int npatterns = 1,
                                const int ncores    = 1) {
  const LMSModel M(modelR);
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);
  return completeLogLikFromModelAghq(M, Z, Gamma, data, colidx, n, npatterns, ncores);
}


// M-step gradient: reverse-mode adjoint of the complete-data log-lik w.r.t.
// theta, accumulated per (row, node) pair via completeLogLikScore (nu=y_i,
// S=0) and accumulateMuSigmaAdjoints -- identical building blocks to the
// shared-grid path (completeGradientReverseFromModel), just re-looped.
// Internal version (see completeLogLikFromModelAghq for why); `setThreads`
// mirrors completeGradientReverseFromModel's flag of the same name -- pass
// false when called from within an already-parallel region (e.g. per-FD-step
// inside the Hessian wrapper) to avoid touching the global OpenMP thread count.
inline arma::vec completeGradientReverseFromModelAghq(
    const LMSModel&  M,
    const arma::cube& Z,
    const arma::mat&  Gamma,
    const std::vector<arma::mat>&  data,
    const std::vector<arma::uvec>& colidx,
    const arma::uvec& n,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric,
    const int npatterns = 1,
    const int ncores    = 1,
    const bool setThreads = true) {
  const arma::uword J    = Z.n_cols;
  const arma::uword k    = Z.n_slices;
  const arma::uword pObs = M.lX.n_rows;
  const std::size_t npar = block.n_elem;

  bool failed = false;

#ifdef _OPENMP
  int oldThreads = omp_get_max_threads();
  if (setThreads) {
    if (ncores <= 0) Rcpp::stop("ncores must be positive");
    omp_set_num_threads(ncores);
  }
  const bool inParallel = omp_in_parallel();
  const int nthreads = inParallel ? 1 : omp_get_max_threads();
#else
  const int nthreads = 1;
#endif

  std::vector<LMSAdjoints> threadAdj;
  threadAdj.reserve(nthreads);
  for (int t = 0; t < nthreads; ++t) threadAdj.emplace_back(M);

  // Per-thread scratch buffers for the padded-to-pObs adjoint inputs, reused
  // across every (row, node) pair instead of allocated fresh each time (that
  // was measured to be a real cost at N*J scale). Entries outside `idx` are
  // never written except via the explicit re-zero below when `idx` changes
  // between patterns, and entries inside `idx` are fully overwritten (not
  // accumulated) on every use, so a single zero per (pattern, thread) is
  // sufficient -- not per (row, node).
  std::vector<arma::vec> threadMuBarFull;
  std::vector<arma::mat> threadSigBarFull;
  threadMuBarFull.reserve(nthreads);
  threadSigBarFull.reserve(nthreads);
  for (int t = 0; t < nthreads; ++t) {
    threadMuBarFull.emplace_back(pObs);
    threadSigBarFull.emplace_back(pObs, pObs);
  }

  arma::uword offset = 0;
  for (int pat = 0; pat < npatterns; ++pat) {
    const arma::uvec& idx  = colidx[pat];
    const arma::mat&  Dat  = data[pat];
    const arma::uword nPat = n[pat];
    const arma::mat Szero  = arma::zeros<arma::mat>(idx.n_elem, idx.n_elem);

    for (int t = 0; t < nthreads; ++t) {
      threadMuBarFull[t].zeros();
      threadSigBarFull[t].zeros();
    }

#pragma omp parallel for default(none) if(!inParallel) schedule(static) \
    shared(M, Dat, idx, Z, Gamma, offset, nPat, J, k, pObs, threadAdj, Szero, \
           threadMuBarFull, threadSigBarFull) \
    reduction(||:failed)
    for (arma::uword i = 0; i < nPat; ++i) {
      const arma::uword rowIdx = offset + i;
      const arma::vec y = Dat.row(i).t();

#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      arma::vec& muBarFull  = threadMuBarFull[tid];
      arma::mat& SigBarFull = threadSigBarFull[tid];

      for (arma::uword j = 0; j < J; ++j) {
        const double tg = Gamma(rowIdx, j);
        if (tg <= DBL_MIN) continue;

        arma::vec zj(k);
        for (arma::uword l = 0; l < k; ++l) zj[l] = Z(rowIdx, j, l);

        const LMSModel::ForwardCache fc = M.muSigmaForward(zj);

        arma::vec muBar_i;
        arma::mat SigBar_i;
        if (!completeLogLikScore(fc.mu.elem(idx), fc.Sigma.submat(idx, idx), y, Szero,
                                 tg, muBar_i, SigBar_i)) {
          failed = true;
          continue;
        }

        muBarFull.elem(idx) = muBar_i;
        SigBarFull.submat(idx, idx) = SigBar_i;

        accumulateMuSigmaAdjointsFromCache(M, fc, muBarFull, SigBarFull, threadAdj[tid]);
      }
    }
    offset += nPat;
  }

  if (failed) {
    arma::vec grad(npar);
    grad.fill(arma::datum::nan);
#ifdef _OPENMP
    if (setThreads) omp_set_num_threads(oldThreads);
#endif
    return grad;
  }

  LMSAdjoints adj(M);
  for (int t = 0; t < nthreads; ++t) adj.add(threadAdj[t]);

  arma::vec grad(npar, arma::fill::zeros);
  for (std::size_t kk = 0; kk < npar; ++kk) {
    const arma::mat& Ablk = lmsAdjointBlock(adj, block[kk]);
    if (Ablk.is_empty()) continue;
    grad[kk] = Ablk(row[kk], col[kk]);
    if (symmetric[kk] && row[kk] != col[kk]) grad[kk] += Ablk(col[kk], row[kk]);
  }

#ifdef _OPENMP
  if (setThreads) omp_set_num_threads(oldThreads);
#endif

  return grad;
}


// [[Rcpp::export]]
arma::vec gradLogLikLmsAghqCpp(const Rcpp::List& modelR,
                               const arma::cube& Z,
                               const arma::mat&  Gamma,
                               const Rcpp::List& dataR,
                               const Rcpp::List& colidxR,
                               const arma::uvec& n,
                               const arma::uvec& block,
                               const arma::uvec& row,
                               const arma::uvec& col,
                               const arma::uvec& symmetric,
                               const int npatterns = 1,
                               const int ncores    = 1) {
  const LMSModel M(modelR);
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);
  return completeGradientReverseFromModelAghq(M, Z, Gamma, data, colidx, n,
                                              block, row, col, symmetric,
                                              npatterns, ncores, true);
}


inline arma::vec completeGradientReverseFromModel(
    const LMSModel&  M,
    const arma::mat& V,
    const std::vector<arma::vec>& TGamma,
    const std::vector<std::vector<arma::vec>>& MeanPatterns,
    const std::vector<std::vector<arma::mat>>& CovPatterns,
    const std::vector<arma::uvec>& colidx,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric,
    const int npatterns = 1,
    const int ncores = 1L,
    const bool setThreads = true) {

  LMSAdjoints adj(M);
  const std::size_t J = V.n_rows;
  const arma::uword pObs = M.lX.n_rows;
  const std::size_t npar = block.n_elem;
  bool failed = false;

#ifdef _OPENMP
  int oldThreads = omp_get_max_threads();
  if (setThreads) {
    if (ncores <= 0)
      Rcpp::stop("ncores must be positive");
    omp_set_num_threads(ncores);
  }
  const bool inParallel = omp_in_parallel();
  const int nthreads = inParallel ? 1 : omp_get_max_threads();
#else
  const int nthreads = 1;
#endif
  std::vector<LMSAdjoints> threadAdj;
  threadAdj.reserve(nthreads);
  for (int t = 0; t < nthreads; ++t)
    threadAdj.emplace_back(M);

#pragma omp parallel for default(none) if(!inParallel) \
  shared(M, V, TGamma, MeanPatterns, CovPatterns, colidx, J, pObs, npatterns, \
         threadAdj, inParallel) reduction(||:failed) schedule(static)
  for (std::size_t j = 0; j < J; j++) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    const arma::vec z = V.row(j).t();
    const LMSModel::ForwardCache fc = M.muSigmaForward(z);

    arma::vec muBarFull(pObs);
    arma::mat SigBarFull(pObs, pObs);
    muBarFull.zeros();
    SigBarFull.zeros();

    for (int i = 0; i < npatterns; i++) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;

      arma::vec muBar_i;
      arma::mat SigBar_i;
      const arma::uvec& idx = colidx[i];

      if (!completeLogLikScore(fc.mu.elem(idx),
                               fc.Sigma.submat(idx, idx),
                               MeanPatterns[j][i],
                               CovPatterns[j][i],
                               tg,
                               muBar_i,
                               SigBar_i)) {
        failed = true;
        continue;
      }

      muBarFull.elem(idx) += muBar_i;
      SigBarFull.submat(idx, idx) += SigBar_i;
    }

    accumulateMuSigmaAdjointsFromCache(M, fc, muBarFull, SigBarFull, threadAdj[tid]);
  }

  if (failed) {
    arma::vec grad(npar);
    grad.fill(arma::datum::nan);
#ifdef _OPENMP
    if (setThreads) omp_set_num_threads(oldThreads);
#endif
    return grad;
  }

  for (int t = 0; t < nthreads; ++t)
    adj.add(threadAdj[t]);

  arma::vec grad(npar, arma::fill::zeros);
  for (std::size_t k = 0; k < npar; ++k) {
    const arma::mat& A = lmsAdjointBlock(adj, block[k]);
    if (A.is_empty()) continue;

    grad[k] = A(row[k], col[k]);
    if (symmetric[k] && row[k] != col[k])
      grad[k] += A(col[k], row[k]);
  }

#ifdef _OPENMP
  if (setThreads) omp_set_num_threads(oldThreads);
#endif

  return grad;
}


// Per-observation complete-data scores at one quadrature node. The node-level
// mean/covariance factorization is shared across all rows in each missing-data
// pattern; only the observation-specific MVN adjoints are recomputed.
// [[Rcpp::export]]
arma::mat completeScoresNodeAnalyticalLmsCpp(
    const Rcpp::List& modelR,
    const Rcpp::List& dataR,
    const arma::vec& z,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric,
    const Rcpp::List& colidxR,
    const arma::uvec& n,
    const int npatterns = 1L,
    const int ncores = 1L) {

  const LMSModel M(modelR);
  const auto data = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);
  // z is shared by every row in this call, so the forward pass (and its
  // adjoint-reverse-pass inputs) only needs computing once, not once per row.
  const LMSModel::ForwardCache fc = M.muSigmaForward(z);
  const arma::vec& mu = fc.mu;
  const arma::mat& Sig = fc.Sigma;
  const arma::uword pObs = M.lX.n_rows;
  const arma::uword npar = block.n_elem;
  const arma::uword N = arma::sum(n);
  arma::mat scores(N, npar, arma::fill::zeros);
  bool failed = false;

  ThreadSetter ts(ncores);
  arma::uword offset = 0;
  for (int pat = 0; pat < npatterns; ++pat) {
    const arma::uvec& idx = colidx[pat];
    const arma::vec muS = mu.elem(idx);
    const arma::mat sigS = Sig.submat(idx, idx);
    arma::mat L;
    if (!sigS.is_finite() || !arma::chol(L, sigS, "lower")) { // avoid warning in arma::chol
      scores.fill(arma::datum::nan);
      return scores;
    }

    const arma::mat Sinv = arma::solve(
      arma::trimatu(L.t()),
      arma::solve(arma::trimatl(L),
                  arma::eye<arma::mat>(sigS.n_rows, sigS.n_cols),
                  arma::solve_opts::fast),
      arma::solve_opts::fast
    );
    const arma::mat centered = data[pat].each_row() - muS.t();
    const arma::mat invDiff = Sinv * centered.t();
    const arma::uword nPat = n[pat];

#pragma omp parallel for default(none) if(ncores > 1) \
    shared(M, fc, idx, Sinv, invDiff, block, row, col, symmetric, scores, \
           offset, nPat, npar, pObs) reduction(||:failed) schedule(static)
    for (arma::uword i = 0; i < nPat; ++i) {
      const arma::vec muBar_i = invDiff.col(i);
      const arma::mat SigmaBar_i = 0.5 *
        (muBar_i * muBar_i.t() - Sinv);
      arma::vec muBar(pObs);
      arma::mat SigmaBar(pObs, pObs);
      muBar.zeros();
      SigmaBar.zeros();
      muBar.elem(idx) = muBar_i;
      SigmaBar.submat(idx, idx) = SigmaBar_i;

      LMSAdjoints adj(M);
      accumulateMuSigmaAdjointsFromCache(M, fc, muBar, SigmaBar, adj);

      for (arma::uword k = 0; k < npar; ++k) {
        const arma::mat& A = lmsAdjointBlock(adj, block[k]);
        if (A.is_empty()) continue;
        double value = A(row[k], col[k]);
        if (symmetric[k] && row[k] != col[k])
          value += A(col[k], row[k]);
        scores(offset + i, k) = value;
        if (!std::isfinite(value)) failed = true;
      }
    }
    offset += nPat;
  }

  if (failed) scores.fill(arma::datum::nan);
  return scores;
}


// Louis' identity building blocks for AGHQ: accumulates, in one pass over
// (row, node) pairs, total_M = sum_i sum_j r_ij * s_ij s_ij' (raw "locations"
// space, npar x npar) and Sbar[i,] = sum_j (r_ij / sqrt(sw_i)) * s_ij, where
// s_ij is the RAW (tgamma=1) complete-data score at that row's own adapted
// node. Mirrors completeScoresNodeAnalyticalLmsCpp's per-observation adjoint
// machinery just above, but looped per-row (own z) instead of per-shared-node
// (one z for every row) -- there's no cross-row sharing to exploit once nodes
// are row-specific, so this accumulates directly instead of returning a
// per-node score matrix for the R side to combine (as the shared-node path
// does). The R side maps total_M/Sbar through the per-group Jacobian into
// theta-space once, mirroring how completeScoresNodeAnalytical() does that
// mapping for the shared-node path today.
// [[Rcpp::export]]
Rcpp::List louisRawScoresAghqCpp(const Rcpp::List& modelR,
                                 const arma::cube& Z,
                                 const arma::mat&  Gamma,
                                 const arma::vec&  sqrtSamplingWeights,
                                 const Rcpp::List& dataR,
                                 const Rcpp::List& colidxR,
                                 const arma::uvec& n,
                                 const arma::uvec& block,
                                 const arma::uvec& row,
                                 const arma::uvec& col,
                                 const arma::uvec& symmetric,
                                 const int npatterns = 1,
                                 const int ncores    = 1) {
  const LMSModel M(modelR);
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);
  const arma::uword J    = Z.n_cols;
  const arma::uword k    = Z.n_slices;
  const arma::uword pObs = M.lX.n_rows;
  const std::size_t npar = block.n_elem;
  const arma::uword N    = arma::sum(n);

  arma::mat Sbar(N, npar, arma::fill::zeros);

  ThreadSetter ts(ncores);
#ifdef _OPENMP
  const int nthreads = omp_get_max_threads();
#else
  const int nthreads = 1;
#endif
  std::vector<arma::mat> threadM;
  threadM.reserve(nthreads);
  for (int t = 0; t < nthreads; ++t)
    threadM.emplace_back(npar, npar, arma::fill::zeros);

  // Per-thread scratch, reused across (row, node) pairs -- see the identical
  // comment in completeGradientReverseFromModelAghq for why a single zero
  // per (pattern, thread) suffices for the pObs-sized buffers. s_ij is fully
  // overwritten (or explicitly zeroed) every use regardless, so it's reused
  // purely to avoid the allocation, not the zeroing.
  std::vector<arma::vec> threadMuBarFull;
  std::vector<arma::mat> threadSigBarFull;
  std::vector<arma::vec> threadSij;
  threadMuBarFull.reserve(nthreads);
  threadSigBarFull.reserve(nthreads);
  threadSij.reserve(nthreads);
  for (int t = 0; t < nthreads; ++t) {
    threadMuBarFull.emplace_back(pObs);
    threadSigBarFull.emplace_back(pObs, pObs);
    threadSij.emplace_back(npar);
  }

  bool failed = false;
  arma::uword offset = 0;

  for (int pat = 0; pat < npatterns; ++pat) {
    const arma::uvec& idx  = colidx[pat];
    const arma::mat&  Dat  = data[pat];
    const arma::uword nPat = n[pat];
    const arma::mat Szero  = arma::zeros<arma::mat>(idx.n_elem, idx.n_elem);

    for (int t = 0; t < nthreads; ++t) {
      threadMuBarFull[t].zeros();
      threadSigBarFull[t].zeros();
    }

#pragma omp parallel for default(none) schedule(static) if(ncores > 1) \
    shared(M, Dat, idx, Z, Gamma, sqrtSamplingWeights, offset, nPat, J, k, \
           pObs, npar, block, row, col, symmetric, Szero, threadM, Sbar, \
           threadMuBarFull, threadSigBarFull, threadSij) \
    reduction(||:failed)
    for (arma::uword i = 0; i < nPat; ++i) {
      const arma::uword rowIdx = offset + i;
      const arma::vec y = Dat.row(i).t();
      const double sqrtSw = sqrtSamplingWeights[rowIdx];

#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      arma::vec& muBarFull  = threadMuBarFull[tid];
      arma::mat& SigBarFull = threadSigBarFull[tid];
      arma::vec& s_ij       = threadSij[tid];

      for (arma::uword j = 0; j < J; ++j) {
        const double r_ij = Gamma(rowIdx, j);
        if (r_ij <= 0.0) continue;

        arma::vec zj(k);
        for (arma::uword l = 0; l < k; ++l) zj[l] = Z(rowIdx, j, l);

        const LMSModel::ForwardCache fc = M.muSigmaForward(zj);

        arma::vec muBar_i;
        arma::mat SigBar_i;
        // raw (unweighted) score: tgamma = 1, S = 0
        if (!completeLogLikScore(fc.mu.elem(idx), fc.Sigma.submat(idx, idx), y, Szero,
                                 1.0, muBar_i, SigBar_i)) {
          failed = true;
          continue;
        }

        muBarFull.elem(idx) = muBar_i;
        SigBarFull.submat(idx, idx) = SigBar_i;

        LMSAdjoints adj(M);
        accumulateMuSigmaAdjointsFromCache(M, fc, muBarFull, SigBarFull, adj);

        s_ij.zeros();
        for (std::size_t kk = 0; kk < npar; ++kk) {
          const arma::mat& Ablk = lmsAdjointBlock(adj, block[kk]);
          if (Ablk.is_empty()) continue;
          double v = Ablk(row[kk], col[kk]);
          if (symmetric[kk] && row[kk] != col[kk]) v += Ablk(col[kk], row[kk]);
          s_ij[kk] = v;
        }

        threadM[tid] += r_ij * (s_ij * s_ij.t());

        if (sqrtSw > 0.0)
          Sbar.row(rowIdx) += (r_ij / sqrtSw) * s_ij.t();
      }
    }
    offset += nPat;
  }

  arma::mat total_M(npar, npar, arma::fill::zeros);
  for (int t = 0; t < nthreads; ++t) total_M += threadM[t];

  if (failed) {
    total_M.fill(arma::datum::nan);
    Sbar.fill(arma::datum::nan);
  }

  return Rcpp::List::create(
    Rcpp::Named("total_M") = total_M,
    Rcpp::Named("Sbar")    = Sbar
  );
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
    LMSModel Mc = M.threadClone();

    // Access the parameter(s) to perturb in the *local* model:
    double& ti   = lmsParam(Mc, block[k], row[k], col[k]);
    double* tj   = nullptr;

    if (symmetric[k] && row[k] != col[k]) {
      tj = &lmsParam(Mc, block[k], col[k], row[k]); // symmetric partner
    }

    // Forward finite difference step
    ti += eps;
    if (tj) *tj += eps;
    Mc.updateCache();

    // Evaluate on the perturbed *local* model
    const double f1 = logLik(Mc);

    // Gradient component
    grad[k] = (f1 - f0) / eps;

    // No need to restore: Mc is thread-local and will be destroyed here.
  }

  return grad;
}


template<class F>
void gradientFDSelected(LMSModel&         M,
                        F&&               logLik,
                        const arma::uvec& block,
                        const arma::uvec& row,
                        const arma::uvec& col,
                        const arma::uvec& symmetric,
                        const arma::uvec& selected,
                        arma::vec&        grad,
                        const double      eps = 1e-6,
                        const int         ncores = 1L) {
  if (selected.is_empty()) return;

  ThreadSetter ts(ncores);
  const double f0 = logLik(M);

  #pragma omp parallel for default(none) \
      shared(M, block, row, col, symmetric, eps, grad, f0, selected) \
      firstprivate(logLik) \
      schedule(static)
  for (arma::uword s = 0; s < selected.n_elem; ++s) {
    const std::size_t k = selected[s];
    LMSModel Mc = M.threadClone();

    double& ti = lmsParam(Mc, block[k], row[k], col[k]);
    double* tj = nullptr;

    if (symmetric[k] && row[k] != col[k]) {
      tj = &lmsParam(Mc, block[k], col[k], row[k]);
    }

    ti += eps;
    if (tj) *tj += eps;
    Mc.updateCache();

    const double f1 = logLik(Mc);
    grad[k] = (f1 - f0) / eps;
  }
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

    const arma::vec& z = V.row(j).t();   // view – no copy
    const std::pair<arma::vec, arma::mat> ms = M.muSigma(z);
    const arma::vec& mu  = ms.first;
    const arma::mat& Sig = ms.second;

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
  const LMSModel model(modelR);

  const arma::mat V      = Rcpp::as<arma::mat>(P["V"]);
  const auto TGamma      = as_vec_of_vec(P["tgamma"]);
  const auto Mean        = as_vec_of_vec_of_vec(P["mean"]);
  const auto Cov         = as_vec_of_vec_of_mat(P["cov"]);
  const auto colidx      = as_vec_of_uvec(colidxR);

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

  const arma::mat V      = Rcpp::as<arma::mat>(P["V"]);
  const auto TGamma      = as_vec_of_vec(P["tgamma"]);
  const auto Mean        = as_vec_of_vec_of_vec(P["mean"]);
  const auto Cov         = as_vec_of_vec_of_mat(P["cov"]);
  const auto colidx      = as_vec_of_uvec(colidxR);

  return completeGradientReverseFromModel(M, V, TGamma, Mean, Cov, colidx,
                                          block, row, col, symmetric,
                                          npatterns, ncores);
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
  const int         N = (int)arma::sum(n);

  ThreadSetter ts(ncores);

#ifdef _OPENMP
  const int nthreads = omp_get_max_threads();
#else
  const int nthreads = 1;
#endif

  // Per-thread density accumulators — avoids atomic writes on the shared vector.
  std::vector<arma::vec> localDens(nthreads, arma::zeros<arma::vec>(N));

#pragma omp parallel for schedule(static) default(none) if(ncores > 1) \
    shared(M, V, w, data, colidx, n, Q, npatterns, localDens)
  for (std::size_t i = 0; i < Q; ++i) {
    if (w[i] <= DBL_MIN) continue;

#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    const arma::vec  z   = V.row(i).t();
    const auto       ms  = M.muSigma(z);
    const arma::vec& mu  = ms.first;
    const arma::mat& Sig = ms.second;

    int offset = 0;
    for (int j = 0; j < npatterns; ++j) {
      const int end = offset + (int)n[j] - 1;
      localDens[tid].subvec(offset, end) +=
        mvnDensSt(data[j], colidx[j], mu, Sig) * w[i];
      offset = end + 1;
    }
  }

  arma::vec density = arma::zeros<arma::vec>(N);
  for (int t = 0; t < nthreads; ++t) density += localDens[t];

  return arma::sum(samplingWeights % arma::log(density));
}


inline arma::vec observedGradientReverseFromModel(
    const LMSModel& M,
    const arma::mat& V,
    const arma::vec& w,
    const arma::vec& samplingWeights,
    const std::vector<arma::mat>& data,
    const std::vector<arma::uvec>& colidx,
    const arma::uvec& n,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric,
    const int npatterns = 1L,
    const int ncores = 1L) {
  // Fisher's identity: the observed-data score is the posterior expectation
  // of the complete-data score.  Recompute the posterior probabilities at M
  // (using the supplied, fixed quadrature), construct the same sufficient
  // statistics as the E-step, and feed them to the analytical reverse-mode
  // complete-data score.
  const std::size_t Q = V.n_rows;
  const arma::uword N = arma::sum(n);
  arma::mat joint(N, Q, arma::fill::zeros);

#pragma omp parallel for schedule(static) default(none) if(ncores > 1) \
  shared(M, V, w, data, colidx, n, joint, Q, npatterns)
  for (std::size_t j = 0; j < Q; ++j) {
    if (w[j] <= DBL_MIN) continue;

    const arma::vec z = V.row(j).t();
    const auto ms = M.muSigma(z);
    int offset = 0;
    for (int pat = 0; pat < npatterns; ++pat) {
      const int end = offset + static_cast<int>(n[pat]) - 1;
      joint(arma::span(offset, end), j) =
        w[j] * mvnDensSt(data[pat], colidx[pat], ms.first, ms.second);
      offset = end + 1;
    }
  }

  const arma::vec marginal = arma::sum(joint, 1);
  if (!marginal.is_finite() || arma::any(marginal <= 0.0)) {
    arma::vec out(block.n_elem);
    out.fill(arma::datum::nan);
    return out;
  }

  arma::mat posterior = joint.each_col() / marginal;
  posterior.each_col() %= samplingWeights;

  std::vector<std::vector<arma::vec>> Mean(
    Q, std::vector<arma::vec>(npatterns));
  std::vector<std::vector<arma::mat>> Cov(
    Q, std::vector<arma::mat>(npatterns));
  std::vector<arma::vec> TGamma(Q, arma::zeros<arma::vec>(npatterns));

#pragma omp parallel for schedule(static) default(none) if(ncores > 1) \
  shared(posterior, data, n, Mean, Cov, TGamma, Q, npatterns)
  for (std::size_t j = 0; j < Q; ++j) {
    int offset = 0;
    for (int pat = 0; pat < npatterns; ++pat) {
      const int end = offset + static_cast<int>(n[pat]) - 1;
      const arma::vec pj = posterior.col(j).subvec(offset, end);
      const double tg = arma::sum(pj);
      TGamma[j][pat] = tg;

      if (tg > DBL_MIN) {
        const arma::vec mean = data[pat].t() * pj / tg;
        const arma::mat centered = data[pat].each_row() - mean.t();
        Mean[j][pat] = mean;
        Cov[j][pat] = centered.t() * (centered.each_col() % pj);
      } else {
        Mean[j][pat] = arma::zeros<arma::vec>(data[pat].n_cols);
        Cov[j][pat] = arma::zeros<arma::mat>(data[pat].n_cols,
                                             data[pat].n_cols);
      }
      offset = end + 1;
    }
  }

  return completeGradientReverseFromModel(
    M, V, TGamma, Mean, Cov, colidx, block, row, col, symmetric,
    npatterns, ncores, false);
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
  ThreadSetter ts(ncores);

  return observedGradientReverseFromModel(
    M, V, w, samplingWeights, data, colidx, n, block, row, col, symmetric,
    npatterns, ncores);
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


inline arma::vec getParams(const LMSModel& M,
                            const arma::uvec& block,
                            const arma::uvec& row,
                            const arma::uvec& col) {
  const std::size_t p = block.n_elem;
  arma::vec pars(p);
  for (std::size_t k = 0; k < p; ++k)
    pars[k] = lmsParam(const_cast<LMSModel&>(M),
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
    double& ti = lmsParam(M, block[k], row[k], col[k]);
    ti = vals[k];

    if (symmetric[k] && row[k] != col[k])
      lmsParam(M, block[k], col[k], row[k]) = vals[k];
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
    LMSModel Mc = M.threadClone();
    setParams(Mc, block, row, col, symmetric, base + disp[k] % incr);
    Mc.updateCache();
    y[k] = fun(Mc);
  }

  // Restore baseline
  setParams(M, block, row, col, symmetric, base);

  // Build design matrix
  const std::size_t q = 1 + 2*p + (p*(p-1))/2;
  arma::mat X(m, q, arma::fill::ones);
  std::size_t colId = 1;
  for (std::size_t j = 0; j < p; ++j, ++colId)
    for (std::size_t k = 0; k < m; ++k)
      X(k, colId) = disp[k][j];
  for (std::size_t j = 0; j < p; ++j, ++colId)
    for (std::size_t k = 0; k < m; ++k)
      X(k, colId) = std::pow(disp[k][j], 2);
  for (std::size_t i = 0; i < p-1; ++i)
    for (std::size_t j = i+1; j < p; ++j, ++colId)
      for (std::size_t k = 0; k < m; ++k)
        X(k, colId) = disp[k][i] * disp[k][j];

  // frac scaling
  arma::vec frac(q, arma::fill::ones);
  for (std::size_t j = 0; j < p; ++j)              frac[1 + j]     = incr[j];
  for (std::size_t j = 0; j < p; ++j)              frac[1 + p + j] = incr[j]*incr[j];
  colId = 1 + 2*p;
  for (std::size_t i = 0; i < p-1; ++i)
    for (std::size_t j = i+1; j < p; ++j, ++colId)
      frac[colId] = incr[i] * incr[j];

  arma::vec coef = arma::solve(X, y) / frac;

  arma::vec grad = coef.subvec(1, p);
  arma::mat Hess(p, p, arma::fill::zeros);
  for (std::size_t j = 0; j < p; ++j)
    Hess(j, j) = 2.0 * coef[1 + p + j];
  colId = 1 + 2*p;
  for (std::size_t i = 0; i < p-1; ++i)
    for (std::size_t j = i+1; j < p; ++j, ++colId) {
      Hess(i,j) = coef[colId];
      Hess(j,i) = coef[colId];
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

  std::vector<std::size_t> idxIp(p), idxIm(p);
  for (std::size_t i=0; i<p; ++i) {
    arma::vec v = arma::zeros<arma::vec>(p);
    v[i]= 1; idxIp[i]=disp.size(); disp.push_back(v);
    v[i]=-1; idxIm[i]=disp.size(); disp.push_back(v);
  }

  std::vector<std::size_t> idxPp(npairs), idxPm(npairs),
                           idxMp(npairs), idxMm(npairs);
  if (p>1) {
    for (std::size_t i=0; i<p-1; ++i)
      for (std::size_t j=i+1; j<p; ++j) {
        std::size_t k = pairIndex(i,j);
        arma::vec v = arma::zeros<arma::vec>(p);
        v[i]= 1; v[j]= 1; idxPp[k]=disp.size(); disp.push_back(v);
        v[i]= 1; v[j]=-1; idxPm[k]=disp.size(); disp.push_back(v);
        v[i]=-1; v[j]= 1; idxMp[k]=disp.size(); disp.push_back(v);
        v[i]=-1; v[j]=-1; idxMm[k]=disp.size(); disp.push_back(v);
      }
  }

  // Evaluate fun (parallel)
  arma::vec y(disp.size());
#pragma omp parallel for default(none) \
  shared(M, disp, block, row, col, symmetric, base, incr, y) \
  firstprivate(fun) schedule(static)
  for (std::size_t k=0; k<disp.size(); ++k) {
    LMSModel Mc = M.threadClone();
    setParams(Mc, block, row, col, symmetric, base + disp[k] % incr);
    Mc.updateCache();
    y[k] = fun(Mc);
  }
  setParams(M, block, row, col, symmetric, base);

  // Assemble gradient/Hessian
  arma::vec grad(p, arma::fill::zeros);
  arma::mat Hess(p, p, arma::fill::zeros);
  const double f0 = y[idx0];

  for (std::size_t i=0; i<p; ++i) {
    double hi = incr[i];
    double fIp = y[idxIp[i]];
    double fIm = y[idxIm[i]];
    grad[i]  = (fIp - fIm) / (2.0*hi);
    Hess(i,i)= (fIp + fIm - 2.0*f0) / (hi*hi);
  }

  if (p>1) {
    for (std::size_t i=0; i<p-1; ++i) {
      double hi = incr[i];
      for (std::size_t j=i+1; j<p; ++j) {
        double hj = incr[j];
        std::size_t k = pairIndex(i,j);
        double fpp=y[idxPp[k]], fpm=y[idxPm[k]],
               fmp=y[idxMp[k]], fmm=y[idxMm[k]];
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

  auto mLs = 1 + 2*p + (p*(p-1))/2;
  auto bytesX = (unsigned long long)mLs * (unsigned long long)mLs *
                (unsigned long long)sizeof(double);

  bool useFullFd = (p >= P_SWITCH) || (bytesX > MEM_LIMIT_BYTES);

  if (!useFullFd)
    return fdHessQuadraticFit(M, std::forward<F>(fun), block, row, col,
        symmetric, base, incr, ncores);
  else
    return fdHessFullFd(M, std::forward<F>(fun), block, row, col,
        symmetric, base, incr, ncores);
}


inline Rcpp::List fdHessFromObservedGradient(
    LMSModel& M,
    const arma::mat& V,
    const arma::vec& w,
    const arma::vec& samplingWeights,
    const std::vector<arma::mat>& data,
    const std::vector<arma::uvec>& colidx,
    const arma::uvec& n,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric,
    const int npatterns,
    const double relStep,
    const double minAbs,
    const int ncores) {

  ThreadSetter ts(ncores);

  const std::size_t p = block.n_elem;
  const arma::vec base = getParams(M, block, row, col);
  const double minScale = (minAbs > 0.0) ? minAbs : 1.0;
  const arma::vec incr =
    arma::max(arma::abs(base), arma::vec(p).fill(minScale)) * relStep;

  const double f0 = observedLogLikFromModel(
    M, V, w, samplingWeights, data, colidx, n, npatterns, ncores);
  const arma::vec grad0 = observedGradientReverseFromModel(
    M, V, w, samplingWeights, data, colidx, n, block, row, col, symmetric,
    npatterns, ncores);

  arma::mat Hess(p, p, arma::fill::zeros);

#pragma omp parallel for default(none) \
  shared(M, V, w, samplingWeights, data, colidx, n, block, row, col, symmetric, \
         npatterns, p, base, incr, grad0, Hess) schedule(static)
  for (std::size_t j = 0; j < p; ++j) {
    LMSModel Mc = M.threadClone();
    arma::vec pars = base;
    pars[j] += incr[j];
    setParams(Mc, block, row, col, symmetric, pars);
    Mc.updateCache();

    const arma::vec gradJ = observedGradientReverseFromModel(
      Mc, V, w, samplingWeights, data, colidx, n, block, row, col,
      symmetric, npatterns, 1L);
    Hess.col(j) = (gradJ - grad0) / incr[j];
  }

  Hess = 0.5 * (Hess + Hess.t());

  return Rcpp::List::create(
    Rcpp::Named("mean")     = f0,
    Rcpp::Named("gradient") = grad0,
    Rcpp::Named("Hessian")  = Hess
  );
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
                               const double      relStep   = 1e-6,
                               const double      minAbs    = 0.0,
                               const int         ncores    = 1L) {
  LMSModel M(modelR);

  const arma::mat V               = Rcpp::as<arma::mat>(P["V"]);
  const arma::vec w               = Rcpp::as<arma::vec>(P["w"]);
  const arma::vec samplingWeights = Rcpp::as<arma::vec>(P["sampling.weights"]);

  const auto colidx = as_vec_of_uvec(colidxR);
  const auto data   = as_vec_of_mat(dataR);

  return fdHessFromObservedGradient(
    M, V, w, samplingWeights, data, colidx, n, block, row, col, symmetric,
    npatterns, relStep, minAbs, ncores);
}


inline Rcpp::List fdHessFromCompleteGradient(
    LMSModel& M,
    const arma::mat& V,
    const std::vector<arma::vec>& TGamma,
    const std::vector<std::vector<arma::vec>>& Mean,
    const std::vector<std::vector<arma::mat>>& Cov,
    const std::vector<arma::uvec>& colidx,
    const arma::uvec& n,
    const arma::uvec& d,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric,
    const int npatterns,
    const double relStep,
    const double minAbs,
    const int ncores) {

  ThreadSetter ts(ncores);

  const std::size_t p = block.n_elem;
  const arma::vec base = getParams(M, block, row, col);
  const double minScale = (minAbs > 0.0) ? minAbs : 1.0;
  const arma::vec incr =
    arma::max(arma::abs(base), arma::vec(p).fill(minScale)) * relStep;

  const double f0 = completeLogLikFromModel(M, V, TGamma, Mean, Cov,
                                            colidx, n, d, npatterns);
  const arma::vec grad0 = completeGradientReverseFromModel(
    M, V, TGamma, Mean, Cov, colidx, block, row, col, symmetric,
    npatterns, ncores, false
  );

  arma::mat Hess(p, p, arma::fill::zeros);

#pragma omp parallel for default(none) \
  shared(M, V, TGamma, Mean, Cov, colidx, n, d, block, row, col, symmetric, \
         npatterns, p, base, incr, grad0, Hess) schedule(static)
  for (std::size_t j = 0; j < p; ++j) {
    LMSModel Mc = M.threadClone();
    arma::vec pars = base;
    pars[j] += incr[j];
    setParams(Mc, block, row, col, symmetric, pars);
    Mc.updateCache();

    const arma::vec gradJ = completeGradientReverseFromModel(
      Mc, V, TGamma, Mean, Cov, colidx, block, row, col, symmetric,
      npatterns, 1L, false
    );

    Hess.col(j) = (gradJ - grad0) / incr[j];
  }

  Hess = 0.5 * (Hess + Hess.t());

  return Rcpp::List::create(
    Rcpp::Named("mean")     = f0,
    Rcpp::Named("gradient") = grad0,
    Rcpp::Named("Hessian")  = Hess
  );
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

  const arma::mat V      = Rcpp::as<arma::mat>(P["V"]);
  const auto TGamma      = as_vec_of_vec(P["tgamma"]);
  const auto Mean        = as_vec_of_vec_of_vec(P["mean"]);
  const auto Cov         = as_vec_of_vec_of_mat(P["cov"]);
  const auto colidx      = as_vec_of_uvec(colidxR);

  return fdHessFromCompleteGradient(M, V, TGamma, Mean, Cov, colidx, n, d,
                                    block, row, col, symmetric, npatterns,
                                    relStep, minAbs, ncores);
}


// AGHQ analogue of fdHessFromCompleteGradient: Hessian of the complete-data
// log-lik via central-difference-free (one-sided) FD of the analytical
// gradient, at FIXED (E-step-time) Z/Gamma -- same recipe, just pointed at
// completeLogLikFromModelAghq/completeGradientReverseFromModelAghq.
inline Rcpp::List fdHessFromCompleteGradientAghq(
    LMSModel&          M,
    const arma::cube&  Z,
    const arma::mat&   Gamma,
    const std::vector<arma::mat>&  data,
    const std::vector<arma::uvec>& colidx,
    const arma::uvec& n,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric,
    const int    npatterns,
    const double relStep,
    const double minAbs,
    const int    ncores) {

  ThreadSetter ts(ncores);

  const std::size_t p = block.n_elem;
  const arma::vec base = getParams(M, block, row, col);
  const double minScale = (minAbs > 0.0) ? minAbs : 1.0;
  const arma::vec incr =
    arma::max(arma::abs(base), arma::vec(p).fill(minScale)) * relStep;

  const double f0 = completeLogLikFromModelAghq(M, Z, Gamma, data, colidx, n,
                                                npatterns, ncores);
  const arma::vec grad0 = completeGradientReverseFromModelAghq(
    M, Z, Gamma, data, colidx, n, block, row, col, symmetric,
    npatterns, ncores, false
  );

  arma::mat Hess(p, p, arma::fill::zeros);

#pragma omp parallel for default(none) \
  shared(M, Z, Gamma, data, colidx, n, block, row, col, symmetric, \
         npatterns, p, base, incr, grad0, Hess) schedule(static)
  for (std::size_t j = 0; j < p; ++j) {
    LMSModel Mc = M.threadClone();
    arma::vec pars = base;
    pars[j] += incr[j];
    setParams(Mc, block, row, col, symmetric, pars);
    Mc.updateCache();

    const arma::vec gradJ = completeGradientReverseFromModelAghq(
      Mc, Z, Gamma, data, colidx, n, block, row, col, symmetric,
      npatterns, 1L, false
    );

    Hess.col(j) = (gradJ - grad0) / incr[j];
  }

  Hess = 0.5 * (Hess + Hess.t());

  return Rcpp::List::create(
    Rcpp::Named("mean")     = f0,
    Rcpp::Named("gradient") = grad0,
    Rcpp::Named("Hessian")  = Hess
  );
}


// [[Rcpp::export]]
Rcpp::List hessCompLogLikLmsAghqCpp(const Rcpp::List& modelR,
                                    const arma::cube& Z,
                                    const arma::mat&  Gamma,
                                    const Rcpp::List& dataR,
                                    const arma::uvec& block,
                                    const arma::uvec& row,
                                    const arma::uvec& col,
                                    const arma::uvec& symmetric,
                                    const Rcpp::List& colidxR,
                                    const arma::uvec& n,
                                    const int         npatterns = 1,
                                    const double      relStep   = 1e-6,
                                    const double      minAbs    = 0.0,
                                    const int         ncores    = 1L) {
  LMSModel M(modelR);
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);

  return fdHessFromCompleteGradientAghq(M, Z, Gamma, data, colidx, n,
                                        block, row, col, symmetric, npatterns,
                                        relStep, minAbs, ncores);
}


// [[Rcpp::export]]
arma::mat densityMatrixLmsCpp(const Rcpp::List& modelR,
                               const arma::mat&  V,
                               const Rcpp::List& dataR,
                               const Rcpp::List& colidxR,
                               const arma::uvec& n,
                               const arma::vec&  samplingWeights,
                               const int npatterns = 1,
                               const int ncores    = 1) {
  const LMSModel M(modelR);
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);

  const std::size_t Q  = V.n_rows;
  const int         N  = (int)arma::sum(n);
  const bool        hasSW = samplingWeights.n_elem == (std::size_t)N;

  arma::mat out(N, Q, arma::fill::zeros);

  ThreadSetter ts(ncores);

#pragma omp parallel for schedule(static) default(none) if(ncores > 1) \
    shared(M, V, data, colidx, n, samplingWeights, out, Q, npatterns, hasSW)
  for (std::size_t i = 0; i < Q; ++i) {
    const arma::vec  z   = V.row(i).t();
    const auto       ms  = M.muSigma(z);
    const arma::vec& mu  = ms.first;
    const arma::mat& Sig = ms.second;

    int offset = 0;
    for (int j = 0; j < npatterns; ++j) {
      const int end = offset + (int)n[j] - 1;
      out(arma::span(offset, end), i) = mvnDensSt(data[j], colidx[j], mu, Sig);
      offset = end + 1;
    }

    if (hasSW)
      out.col(i) = arma::exp(arma::log(out.col(i)) % samplingWeights);
  }

  return out;
}


// [[Rcpp::export]]
Rcpp::List estepSuffStatLmsCpp(const arma::mat&  P,
                                const Rcpp::List& dataR,
                                const arma::uvec& n,
                                const int npatterns = 1,
                                const int ncores    = 1) {
  const auto data   = as_vec_of_mat(dataR);
  const std::size_t Q = P.n_cols;

  // Thread-safe intermediates: R memory allocation (Rcpp::List) is not safe
  // inside a parallel region, so accumulate into std::vector first.
  std::vector<std::vector<arma::vec>> allMeans(Q, std::vector<arma::vec>(npatterns));
  std::vector<std::vector<arma::mat>> allCovs (Q, std::vector<arma::mat>(npatterns));
  std::vector<arma::vec>              allTg   (Q, arma::zeros<arma::vec>(npatterns));

  ThreadSetter ts(ncores);

#pragma omp parallel for schedule(static) default(none) if(ncores > 1) \
    shared(P, data, n, Q, npatterns, allMeans, allCovs, allTg)
  for (std::size_t i = 0; i < Q; ++i) {
    const arma::vec p = P.col(i);
    int offset = 0;

    for (int j = 0; j < npatterns; ++j) {
      const int        end = offset + (int)n[j] - 1;
      const arma::vec   pj = p.subvec(offset, end);
      const arma::mat&  Dj = data[j];

      const double tg = arma::sum(pj);
      arma::vec    wm = Dj.t() * pj / tg;
      arma::mat    X  = Dj.each_row() - wm.t();
      arma::mat    cv = X.t() * (X.each_col() % pj);

      allMeans[i][j] = wm;
      allCovs [i][j] = cv;
      allTg   [i][j] = tg;

      offset = end + 1;
    }
  }

  // Convert to Rcpp::List — R memory allocation must be single-threaded.
  Rcpp::List wMeans(Q), wCovs(Q), tGamma(Q);
  for (std::size_t i = 0; i < Q; ++i) {
    Rcpp::List wMeans_i(npatterns), wCovs_i(npatterns);
    arma::vec  tGamma_i(npatterns);

    for (int j = 0; j < npatterns; ++j) {
      wMeans_i[j]  = allMeans[i][j];
      wCovs_i[j]   = allCovs [i][j];
      tGamma_i[j]  = allTg   [i][j];
    }

    wMeans[i] = wMeans_i;
    wCovs[i]  = wCovs_i;
    tGamma[i] = tGamma_i;
  }

  return Rcpp::List::create(
    Rcpp::Named("mean")   = wMeans,
    Rcpp::Named("cov")    = wCovs,
    Rcpp::Named("tgamma") = tGamma
  );
}
