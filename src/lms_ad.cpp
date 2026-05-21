// Reverse-mode AD gradient and central-FD-of-gradient Hessian for LMS
// log-likelihoods via Stan Math.
// Stan headers MUST come first to register Eigen plugins before any
// other header pulls in Eigen.

#include <stan/math/rev.hpp>  // must precede any other Eigen include
#include <RcppEigen.h>
#include "thread_setter.h"

// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

// Suppress OpenBLAS internal threading inside OMP regions.
// Declared weak so the linker resolves it to nullptr on non-OpenBLAS backends
// (reference BLAS, MKL, Accelerate) — the call site guards against null.
#if defined(__GNUC__) || defined(__clang__)
extern "C" { __attribute__((weak)) void openblas_set_num_threads(int n); }
#else
static inline void openblas_set_num_threads(int) {}
#endif

// ─── type aliases ─────────────────────────────────────────────────────────────
template<typename T>
using EMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using EVec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using EMatd = Eigen::MatrixXd;
using EVecd = Eigen::VectorXd;

// ─── Kronecker kron(A_double, B_T) ───────────────────────────────────────────
template<typename T>
EMat<T> kron_d_T(const EMatd& A, const EMat<T>& B) {
  int ma = A.rows(), na = A.cols(), mb = B.rows(), nb = B.cols();
  EMat<T> R(ma * mb, na * nb);
  R.setZero();
  for (int i = 0; i < ma; ++i)
    for (int j = 0; j < na; ++j)
      if (A(i, j) != 0.0)
        R.block(i * mb, j * nb, mb, nb) = A(i, j) * B;
  return R;
}

// ─── Weighted MVN log-likelihood from sufficient statistics ───────────────────
// Computes: -0.5 * [tg*(d*log2pi + log|Sig|) + trace(Sig^{-1}*S) + tg*(nu-mu)'Sig^{-1}(nu-mu)]
// nu, S are double constants (from E-step); mu, Sig are T (model parameters).
template<typename T>
T total_dmvn_weighted_ad(const EVec<T>& mu, const EMat<T>& Sig,
                          const EVecd& nu,   const EMatd& S,
                          double tg, int d) {
  using stan::math::cholesky_decompose;
  using stan::math::mdivide_left_tri_low;
  using stan::math::mdivide_left_spd;

  EMat<T> L = cholesky_decompose(Sig);

  // log|Sig| = 2 * sum(log(diag(L)))
  T log_det = 2.0 * L.diagonal().array().log().sum();

  // Mahalanobis: tg * ||L^{-1}(nu - mu)||^2
  EVec<T> diff(d);
  for (int i = 0; i < d; ++i) diff(i) = nu(i) - mu(i);
  EVec<T> y = mdivide_left_tri_low(L, diff);
  T mahal   = tg * y.squaredNorm();

  // trace(Sig^{-1} S)
  T trace_term = mdivide_left_spd(Sig, S).trace();

  static constexpr double log2pi = 1.8378770664093455;
  return -0.5 * (tg * (d * log2pi + log_det) + trace_term + mahal);
}

// ─── Subset helpers ──────────────────────────────────────────────────────────
template<typename T>
EVec<T> subset_vec(const EVec<T>& v, const std::vector<int>& idx) {
  EVec<T> out(idx.size());
  for (int k = 0; k < (int)idx.size(); ++k) out(k) = v(idx[k]);
  return out;
}

template<typename T>
EMat<T> subset_sym_mat(const EMat<T>& M, const std::vector<int>& idx) {
  int n = (int)idx.size();
  EMat<T> out(n, n);
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < n; ++b)
      out(a, b) = M(idx[a], idx[b]);
  return out;
}

// ─── R-list parsers (arma-free) ───────────────────────────────────────────────
static std::vector<std::vector<EVecd>>
parse_list_of_list_of_vec(const Rcpp::List& L) {
  int J = L.size();
  std::vector<std::vector<EVecd>> out(J);
  for (int j = 0; j < J; ++j) {
    Rcpp::List inner = L[j];
    int I = inner.size();
    out[j].resize(I);
    for (int i = 0; i < I; ++i)
      out[j][i] = Rcpp::as<EVecd>(inner[i]);
  }
  return out;
}

static std::vector<std::vector<EMatd>>
parse_list_of_list_of_mat(const Rcpp::List& L) {
  int J = L.size();
  std::vector<std::vector<EMatd>> out(J);
  for (int j = 0; j < J; ++j) {
    Rcpp::List inner = L[j];
    int I = inner.size();
    out[j].resize(I);
    for (int i = 0; i < I; ++i)
      out[j][i] = Rcpp::as<EMatd>(inner[i]);
  }
  return out;
}

static std::vector<EVecd>
parse_list_of_vec(const Rcpp::List& L) {
  int J = L.size();
  std::vector<EVecd> out(J);
  for (int j = 0; j < J; ++j)
    out[j] = Rcpp::as<EVecd>(L[j]);
  return out;
}

static std::vector<std::vector<int>>
parse_list_of_ivec(const Rcpp::List& L) {
  int J = L.size();
  std::vector<std::vector<int>> out(J);
  for (int j = 0; j < J; ++j) {
    Rcpp::IntegerVector v = L[j];
    out[j].assign(v.begin(), v.end());
  }
  return out;
}

// ─── Functor: complete log-likelihood in T arithmetic ────────────────────────
struct CompLogLikAD {
  // Base model matrices (all double — loaded once from R)
  EMatd A0, Oxx0, Oex0, Ie0, lY0, lX0, W0, Tmat0, tY0, tX0, Gx0, Ge0;
  EMatd a0, beta00, Psi0, d0, e0;
  int k0, numXis0;
  bool hasComposites0;

  // Free-parameter specification (blk, r, c, sym: 0-based)
  std::vector<int> blk, r, c, sym;

  // Sufficient statistics — indexed [j][i] for node j, pattern i
  int Q, npatterns;
  EMatd V;                                       // Q × k quadrature nodes
  std::vector<std::vector<EVecd>> nu;            // weighted means
  std::vector<std::vector<EMatd>> S;             // weighted scatter
  std::vector<EVecd> tg;                         // total weights (Q × npatterns)
  std::vector<std::vector<int>> colidx;          // pattern column indices (0-based)
  std::vector<int> dims;                         // observed dimension per pattern

  // Called by stan::math::gradient (T=var)
  template<typename T>
  T operator()(const EVec<T>& theta) const {
    // ── Build T-typed matrices from base doubles, overriding free params ──
    EMat<T> A     = A0.cast<T>();
    EMat<T> Oxx   = Oxx0.cast<T>();
    EMat<T> Oex   = Oex0.cast<T>();
    EMat<T> lX    = lX0.cast<T>();
    EMat<T> lY    = lY0.cast<T>();
    EMat<T> W     = W0.cast<T>();
    EMat<T> Tmat  = Tmat0.cast<T>();
    EMat<T> tX    = tX0.cast<T>();
    EMat<T> tY    = tY0.cast<T>();
    EMat<T> Gx    = Gx0.cast<T>();
    EMat<T> Ge    = Ge0.cast<T>();
    EMat<T> a     = a0.cast<T>();
    EMat<T> beta0 = beta00.cast<T>();
    EMat<T> Psi   = Psi0.cast<T>();
    EMat<T> dMat  = d0.cast<T>();
    EMat<T> e     = e0.cast<T>();

    for (int ki = 0; ki < (int)theta.size(); ++ki) {
      EMat<T>* mat = nullptr;
      switch (blk[ki]) {
        case 0:  mat = &lX;    break;
        case 1:  mat = &lY;    break;
        case 2:  mat = &tX;    break;
        case 3:  mat = &tY;    break;
        case 4:  mat = &dMat;  break;
        case 5:  mat = &e;     break;
        case 6:  mat = &A;     break;
        case 7:  mat = &Psi;   break;
        case 8:  mat = &a;     break;
        case 9:  mat = &beta0; break;
        case 10: mat = &Gx;    break;
        case 11: mat = &Ge;    break;
        case 12: mat = &Oxx;   break;
        case 13: mat = &Oex;   break;
        case 14: mat = &W;     break;
        case 15: mat = &Tmat;  break;
        default: break;
      }
      if (mat) {
        (*mat)(r[ki], c[ki]) = theta[ki];
        if (sym[ki] && r[ki] != c[ki])
          (*mat)(c[ki], r[ki]) = theta[ki];
      }
    }

    // ── z-independent cached quantities ──────────────────────────────────
    EMatd Oi_d = EMatd::Identity(numXis0, numXis0);
    for (int i = 0; i < k0; ++i) Oi_d(i, i) = 0.0;
    EMat<T> Oi    = Oi_d.cast<T>();
    EMat<T> GxA   = Gx * A;
    EMat<T> AOi   = A * Oi;
    EMat<T> AOiAt = AOi * A.transpose();

    // ── Composite-adjusted loading matrix and residual (z-independent) ───
    // Mirrors muSigma() in equations_lms.cpp:
    //   lXc = lX + T * W * pinv(W^T * T * W)
    //   dc  = d + T - lXc * (W^T * T * W) * lXc^T
    // generalized_inverse (SVD-based) matches arma::pinv and stays valid
    // when WtTW is rank-deficient at perturbed parameter values.
    EMat<T> lX_eff, d_eff;
    if (hasComposites0) {
      EMat<T> WtTW = W.transpose() * Tmat * W;
      lX_eff = lX + Tmat * W * stan::math::generalized_inverse(WtTW);
      d_eff  = dMat + Tmat - lX_eff * WtTW * lX_eff.transpose();
    } else {
      lX_eff = lX;
      d_eff  = dMat;
    }

    // ── Accumulate complete log-likelihood ────────────────────────────────
    T ll = T(0.0);

    for (int j = 0; j < Q; ++j) {
      EVecd z_d = V.row(j).transpose();
      EVecd zVec_d = EVecd::Zero(numXis0);
      if (k0 > 0) zVec_d.head(k0) = z_d;

      EMat<T> v    = beta0 + A * zVec_d.cast<T>();
      EMat<T> kronZ = kron_d_T(Ie0, v);

      EMat<T> B    = Ie0.cast<T>() - Ge - kronZ.transpose() * Oex;
      EMat<T> Binv = stan::math::inverse(B);

      EMat<T>& muXi = v;
      EMat<T> muEta = Binv * (a + Gx * muXi + kronZ.transpose() * Oxx * muXi);
      EMat<T> Eta   = Binv * (GxA + kronZ.transpose() * Oxx * A);

      int nxi = muXi.rows(), neta = muEta.rows();
      EMat<T> varEta   = Eta * Oi * Eta.transpose() + Binv * Psi * Binv.transpose();
      EMat<T> covXiEta = AOi * Eta.transpose();

      EMat<T> vcov(nxi + neta, nxi + neta);
      vcov.topLeftCorner(nxi, nxi)       = AOiAt;
      vcov.topRightCorner(nxi, neta)     = covXiEta;
      vcov.bottomLeftCorner(neta, nxi)   = covXiEta.transpose();
      vcov.bottomRightCorner(neta, neta) = varEta;

      EVec<T> xieta(nxi + neta);
      xieta.head(nxi)  = muXi;
      xieta.tail(neta) = muEta;

      EVec<T> mu_j  = tX + lX_eff * xieta;
      EMat<T> Sig_j = lX_eff * vcov * lX_eff.transpose() + d_eff;

      for (int i = 0; i < npatterns; ++i) {
        if (tg[j](i) <= 0.0) continue;
        EVec<T> mu_i  = subset_vec   (mu_j,  colidx[i]);
        EMat<T> Sig_i = subset_sym_mat(Sig_j, colidx[i]);
        ll += total_dmvn_weighted_ad(mu_i, Sig_i,
                                     nu[j][i], S[j][i],
                                     tg[j](i), dims[i]);
      }
    }

    return ll;
  }
};


// ─── Shared setup: populate functor + build theta0 from R args ───────────────
static std::pair<CompLogLikAD, EVecd>
buildCompLogLikAD(const Rcpp::List&    modelR,
                  const Rcpp::List&    P,
                  const Rcpp::IntegerVector& block,
                  const Rcpp::IntegerVector& row,
                  const Rcpp::IntegerVector& col,
                  const Rcpp::IntegerVector& symmetric,
                  const Rcpp::List&    colidxR,
                  const Rcpp::IntegerVector& n,
                  const Rcpp::IntegerVector& d,
                  const int            npatterns) {
  CompLogLikAD f;

  const Rcpp::List matrices = modelR["matrices"];
  const Rcpp::List info     = modelR["info"];
  const Rcpp::List quad     = modelR["quad"];

  f.A0     = Rcpp::as<EMatd>(matrices["A"]);
  f.Oxx0   = Rcpp::as<EMatd>(matrices["omegaXiXi"]);
  f.Oex0   = Rcpp::as<EMatd>(matrices["omegaEtaXi"]);
  f.Ie0    = Rcpp::as<EMatd>(matrices["Ieta"]);
  f.lY0    = Rcpp::as<EMatd>(matrices["lambdaY"]);
  f.lX0    = Rcpp::as<EMatd>(matrices["lambdaX"]);
  f.tY0    = Rcpp::as<EMatd>(matrices["tauY"]);
  f.tX0    = Rcpp::as<EMatd>(matrices["tauX"]);
  f.W0     = Rcpp::as<EMatd>(matrices["W"]);
  f.Tmat0  = Rcpp::as<EMatd>(matrices["T"]);
  f.Gx0    = Rcpp::as<EMatd>(matrices["gammaXi"]);
  f.Ge0    = Rcpp::as<EMatd>(matrices["gammaEta"]);
  f.a0     = Rcpp::as<EMatd>(matrices["alpha"]);
  f.beta00 = Rcpp::as<EMatd>(matrices["beta0"]);
  f.Psi0   = Rcpp::as<EMatd>(matrices["psi"]);
  f.d0     = Rcpp::as<EMatd>(matrices["thetaDelta"]);
  f.e0     = Rcpp::as<EMatd>(matrices["thetaEpsilon"]);

  f.k0             = Rcpp::as<int>(quad["k"]);
  f.numXis0        = Rcpp::as<int>(info["numXis"]);
  f.hasComposites0 = Rcpp::as<bool>(info["hasComposites"]);

  const int p = block.size();
  f.blk.resize(p); f.r.resize(p); f.c.resize(p); f.sym.resize(p);
  for (int ki = 0; ki < p; ++ki) {
    f.blk[ki] = block[ki];
    f.r  [ki] = row  [ki];
    f.c  [ki] = col  [ki];
    f.sym[ki] = symmetric[ki];
  }

  f.V = Rcpp::as<EMatd>(P["V"]);
  f.Q = f.V.rows();
  f.npatterns = npatterns;

  f.nu = parse_list_of_list_of_vec(P["mean"]);
  f.S  = parse_list_of_list_of_mat(P["cov"]);

  const auto TGamma = parse_list_of_vec(P["tgamma"]);
  f.tg.resize(f.Q);
  for (int j = 0; j < f.Q; ++j) f.tg[j] = TGamma[j];

  f.colidx = parse_list_of_ivec(colidxR);
  f.dims.resize(npatterns);
  for (int i = 0; i < npatterns; ++i) f.dims[i] = d[i];

  // Extract current parameter values for theta0
  auto get_elem = [&](int bi, int ri, int ci) -> double {
    switch (bi) {
      case 0:  return f.lX0   (ri, ci);
      case 1:  return f.lY0   (ri, ci);
      case 2:  return f.tX0   (ri, ci);
      case 3:  return f.tY0   (ri, ci);
      case 4:  return f.d0    (ri, ci);
      case 5:  return f.e0    (ri, ci);
      case 6:  return f.A0    (ri, ci);
      case 7:  return f.Psi0  (ri, ci);
      case 8:  return f.a0    (ri, ci);
      case 9:  return f.beta00(ri, ci);
      case 10: return f.Gx0   (ri, ci);
      case 11: return f.Ge0   (ri, ci);
      case 12: return f.Oxx0  (ri, ci);
      case 13: return f.Oex0  (ri, ci);
      case 14: return f.W0    (ri, ci);
      case 15: return f.Tmat0 (ri, ci);
      default: return 0.0;
    }
  };

  EVecd theta0(p);
  for (int ki = 0; ki < p; ++ki)
    theta0[ki] = get_elem(f.blk[ki], f.r[ki], f.c[ki]);

  return {std::move(f), std::move(theta0)};
}


// ─── AD gradient of complete log-likelihood ──────────────────────────────────
// [[Rcpp::export]]
EVecd gradCompLogLikAdLmsCpp(const Rcpp::List& modelR,
                              const Rcpp::List& P,
                              const Rcpp::IntegerVector& block,
                              const Rcpp::IntegerVector& row,
                              const Rcpp::IntegerVector& col,
                              const Rcpp::IntegerVector& symmetric,
                              const Rcpp::List& colidxR,
                              const Rcpp::IntegerVector& n,
                              const Rcpp::IntegerVector& d,
                              const int npatterns = 1) {
  auto [f, theta0] = buildCompLogLikAD(modelR, P, block, row, col, symmetric,
                                        colidxR, n, d, npatterns);
  double ll_val;
  EVecd  grad_val;
  stan::math::gradient(f, theta0, ll_val, grad_val);
  return grad_val;
}


// ─── Hessian via central-FD of AD gradient ────────────────────────────────────
// For each direction i, two reverse-mode AD gradient evaluations at
//   theta0 ± eps*e_i  give an entire column of the Hessian:
//   H[:,i] = (grad(theta0 + eps*e_i) − grad(theta0 − eps*e_i)) / (2*eps)
// Cost: 1 baseline gradient + 2p gradient evals (each ≈ cost of one function
// eval), totalling O(p) — versus O(p²) for FD-of-function approaches.
// Columns are independent: trivially parallelised over OMP threads.
// Returns List(mean, gradient, Hessian) matching hessCompLogLikLmsCpp format.
// [[Rcpp::export]]
Rcpp::List hessCompLogLikAdLmsCpp(const Rcpp::List& modelR,
                                   const Rcpp::List& P,
                                   const Rcpp::IntegerVector& block,
                                   const Rcpp::IntegerVector& row,
                                   const Rcpp::IntegerVector& col,
                                   const Rcpp::IntegerVector& symmetric,
                                   const Rcpp::List& colidxR,
                                   const Rcpp::IntegerVector& n,
                                   const Rcpp::IntegerVector& d,
                                   const int npatterns = 1,
                                   const int ncores    = 1) {
  auto [f, theta0] = buildCompLogLikAD(modelR, P, block, row, col, symmetric,
                                        colidxR, n, d, npatterns);
  const int p = theta0.size();

  // Baseline: ll and gradient at theta0 (single thread, before OMP region)
  double ll_val;
  EVecd  grad0(p);
  stan::math::gradient(f, theta0, ll_val, grad0);

  // Optimal step for central FD of an exact gradient: (eps_mach)^(1/3) ~ 6e-6
  const double eps = std::cbrt(std::numeric_limits<double>::epsilon());

  EMatd H(p, p);
  ThreadSetter ts(ncores);

  // With STAN_THREADS, each OMP thread needs its own ChainableStack initialised
  // before any var operations.  stan::math::gradient uses nested_rev_autodiff
  // internally, so it is safe once the stack exists.
#pragma omp parallel
  {
    stan::math::ChainableStack thread_stack;
    if (openblas_set_num_threads) openblas_set_num_threads(1);

#pragma omp for schedule(static)
    for (int i = 0; i < p; ++i) {
      EVecd tp = theta0, tm = theta0;
      tp(i) += eps;
      tm(i) -= eps;

      double ll_p, ll_m;
      EVecd  gp(p), gm(p);
      stan::math::gradient(f, tp, ll_p, gp);
      stan::math::gradient(f, tm, ll_m, gm);

      H.col(i) = (gp - gm) / (2.0 * eps);
    }
  }

  H = 0.5 * (H + H.transpose());

  return Rcpp::List::create(
    Rcpp::Named("mean")     = ll_val,
    Rcpp::Named("gradient") = grad0,
    Rcpp::Named("Hessian")  = H
  );
}
