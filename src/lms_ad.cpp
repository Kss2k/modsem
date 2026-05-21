// Reverse-mode AD gradient and forward-over-reverse Hessian for LMS
// log-likelihoods via Stan Math.
// Stan headers MUST come first to register Eigen plugins before any
// other header pulls in Eigen.
// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <stan/math/mix.hpp>   // includes rev + fwd — needed for hessian()
#include <RcppArmadillo.h>
#include "utils.h"
#include "mvnorm.h"             // ThreadSetter + OpenMP guards

// Forward declaration: suppress OpenBLAS internal threading inside OMP regions.
// The symbol is present in libopenblas; falls back to a no-op if absent.
extern "C" { void openblas_set_num_threads(int); }

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
EVec<T> subset_vec(const EVec<T>& v, const arma::uvec& idx) {
  EVec<T> out(idx.n_elem);
  for (arma::uword k = 0; k < idx.n_elem; ++k) out(k) = v(idx[k]);
  return out;
}

template<typename T>
EMat<T> subset_sym_mat(const EMat<T>& M, const arma::uvec& idx) {
  int n = (int)idx.n_elem;
  EMat<T> out(n, n);
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < n; ++b)
      out(a, b) = M(idx[a], idx[b]);
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
  EMatd V;                                     // Q × k quadrature nodes
  std::vector<std::vector<EVecd>> nu;          // weighted means
  std::vector<std::vector<EMatd>> S;           // weighted scatter
  std::vector<std::vector<double>> tg;         // total weights
  std::vector<arma::uvec> colidx;              // pattern column indices (arma, 0-based)
  std::vector<int> dims;                       // observed dimension per pattern

  // Called by stan::math::gradient (T=var) and stan::math::hessian (T=fvar<var>)
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

      EVec<T> mu_j  = tX + lX * xieta;
      EMat<T> Sig_j = lX * vcov * lX.transpose() + dMat;

      for (int i = 0; i < npatterns; ++i) {
        if (tg[j][i] <= 0.0) continue;
        EVec<T> mu_i  = subset_vec   (mu_j,  colidx[i]);
        EMat<T> Sig_i = subset_sym_mat(Sig_j, colidx[i]);
        ll += total_dmvn_weighted_ad(mu_i, Sig_i,
                                     nu[j][i], S[j][i],
                                     tg[j][i], dims[i]);
      }
    }

    return ll;
  }
};


// ─── Helpers to convert arma ↔ Eigen ─────────────────────────────────────────
static EMatd arma2eigen(const arma::mat& m) {
  return Eigen::Map<const EMatd>(m.memptr(), m.n_rows, m.n_cols);
}
static EVecd arma2eigenv(const arma::vec& v) {
  return Eigen::Map<const EVecd>(v.memptr(), v.n_elem);
}


// ─── Shared setup: populate functor + build theta0 from R args ───────────────
static std::pair<CompLogLikAD, EVecd>
buildCompLogLikAD(const Rcpp::List& modelR,
                  const Rcpp::List& P,
                  const arma::uvec& block,
                  const arma::uvec& row,
                  const arma::uvec& col,
                  const arma::uvec& symmetric,
                  const Rcpp::List& colidxR,
                  const arma::uvec& n,
                  const arma::uvec& d,
                  const int         npatterns) {
  CompLogLikAD f;

  const Rcpp::List matrices = modelR["matrices"];
  const Rcpp::List info     = modelR["info"];
  const Rcpp::List quad     = modelR["quad"];

  f.A0     = arma2eigen(Rcpp::as<arma::mat>(matrices["A"]));
  f.Oxx0   = arma2eigen(Rcpp::as<arma::mat>(matrices["omegaXiXi"]));
  f.Oex0   = arma2eigen(Rcpp::as<arma::mat>(matrices["omegaEtaXi"]));
  f.Ie0    = arma2eigen(Rcpp::as<arma::mat>(matrices["Ieta"]));
  f.lY0    = arma2eigen(Rcpp::as<arma::mat>(matrices["lambdaY"]));
  f.lX0    = arma2eigen(Rcpp::as<arma::mat>(matrices["lambdaX"]));
  f.tY0    = arma2eigen(Rcpp::as<arma::mat>(matrices["tauY"]));
  f.tX0    = arma2eigen(Rcpp::as<arma::mat>(matrices["tauX"]));
  f.W0     = arma2eigen(Rcpp::as<arma::mat>(matrices["W"]));
  f.Tmat0  = arma2eigen(Rcpp::as<arma::mat>(matrices["T"]));
  f.Gx0    = arma2eigen(Rcpp::as<arma::mat>(matrices["gammaXi"]));
  f.Ge0    = arma2eigen(Rcpp::as<arma::mat>(matrices["gammaEta"]));
  f.a0     = arma2eigen(Rcpp::as<arma::mat>(matrices["alpha"]));
  f.beta00 = arma2eigen(Rcpp::as<arma::mat>(matrices["beta0"]));
  f.Psi0   = arma2eigen(Rcpp::as<arma::mat>(matrices["psi"]));
  f.d0     = arma2eigen(Rcpp::as<arma::mat>(matrices["thetaDelta"]));
  f.e0     = arma2eigen(Rcpp::as<arma::mat>(matrices["thetaEpsilon"]));

  f.k0             = Rcpp::as<int>(quad["k"]);
  f.numXis0        = Rcpp::as<int>(info["numXis"]);
  f.hasComposites0 = Rcpp::as<bool>(info["hasComposites"]);

  if (f.hasComposites0)
    Rcpp::stop("AD path: composite variables not yet supported");

  const int p = (int)block.n_elem;
  f.blk.resize(p); f.r.resize(p); f.c.resize(p); f.sym.resize(p);
  for (int ki = 0; ki < p; ++ki) {
    f.blk[ki] = (int)block[ki];
    f.r  [ki] = (int)row  [ki];
    f.c  [ki] = (int)col  [ki];
    f.sym[ki] = (int)symmetric[ki];
  }

  const arma::mat Vmat = Rcpp::as<arma::mat>(P["V"]);
  f.V = arma2eigen(Vmat);
  f.Q = f.V.rows();
  f.npatterns = npatterns;

  const auto Mean   = as_vec_of_vec_of_vec(P["mean"]);
  const auto Cov    = as_vec_of_vec_of_mat(P["cov"]);
  const auto TGamma = as_vec_of_vec(P["tgamma"]);

  f.nu.resize(f.Q); f.S.resize(f.Q); f.tg.resize(f.Q);
  for (int j = 0; j < f.Q; ++j) {
    f.nu[j].resize(npatterns);
    f.S [j].resize(npatterns);
    f.tg[j].resize(npatterns);
    for (int i = 0; i < npatterns; ++i) {
      f.nu[j][i] = arma2eigenv(Mean[j][i]);
      f.S [j][i] = arma2eigen (Cov [j][i]);
      f.tg[j][i] = TGamma[j][i];
    }
  }

  f.colidx = as_vec_of_uvec(colidxR);
  f.dims.resize(npatterns);
  for (int i = 0; i < npatterns; ++i) f.dims[i] = (int)d[i];

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
arma::vec gradCompLogLikAdLmsCpp(const Rcpp::List& modelR,
                                  const Rcpp::List& P,
                                  const arma::uvec& block,
                                  const arma::uvec& row,
                                  const arma::uvec& col,
                                  const arma::uvec& symmetric,
                                  const Rcpp::List& colidxR,
                                  const arma::uvec& n,
                                  const arma::uvec& d,
                                  const int npatterns = 1) {
  auto [f, theta0] = buildCompLogLikAD(modelR, P, block, row, col, symmetric,
                                        colidxR, n, d, npatterns);
  double ll_val;
  EVecd  grad_val;
  stan::math::gradient(f, theta0, ll_val, grad_val);

  return arma::vec(grad_val.data(), grad_val.size());
}


// ─── AD Hessian of complete log-likelihood ────────────────────────────────────
// Each of the p forward-over-reverse sweeps is independent: they run in
// parallel via OpenMP, each with its own thread-local Stan Math autodiff stack
// (nested_rev_autodiff creates a fresh nested scope per thread).
// Returns List(mean, gradient, Hessian) matching hessCompLogLikLmsCpp format.
// [[Rcpp::export]]
Rcpp::List hessCompLogLikAdLmsCpp(const Rcpp::List& modelR,
                                   const Rcpp::List& P,
                                   const arma::uvec& block,
                                   const arma::uvec& row,
                                   const arma::uvec& col,
                                   const arma::uvec& symmetric,
                                   const Rcpp::List& colidxR,
                                   const arma::uvec& n,
                                   const arma::uvec& d,
                                   const int npatterns = 1,
                                   const int ncores    = 1) {
  auto [f, theta0] = buildCompLogLikAD(modelR, P, block, row, col, symmetric,
                                        colidxR, n, d, npatterns);
  const int p = theta0.size();

  using FVV = stan::math::fvar<stan::math::var>;

  EMatd  H(p, p);
  EVecd  grad_val(p);
  double ll_val      = 0.0;
  bool   ll_captured = false;

  ThreadSetter ts(ncores);

  // With STAN_THREADS, the autodiff stack is TLS (__thread pointer).
  // Each new OMP thread starts with instance_=nullptr and must instantiate
  // a ChainableStack before touching any var/fvar<var> operations.
  // nested_rev_autodiff is then safe: it pushes/pops on the thread's own stack.
  //
  // OpenBLAS (and some other BLAS) will serialize BLAS calls made inside an
  // OMP parallel region unless its own thread count is 1.  Temporarily set it
  // to 1 so that each OMP thread's matrix ops run single-threaded.
#pragma omp parallel
  {
    stan::math::ChainableStack thread_stack;  // initialize TLS for this OMP thread

    // Suppress BLAS internal threading while inside OMP region
    openblas_set_num_threads(1);

#pragma omp for schedule(static)
  for (int i = 0; i < p; ++i) {
    stan::math::nested_rev_autodiff nested;

    EVec<FVV> theta_fvd(p);
    for (int j = 0; j < p; ++j)
      theta_fvd(j) = FVV(theta0(j), j == i ? 1.0 : 0.0);

    FVV result = f(theta_fvd);

    // gradient[i] = d(f)/d(theta_i) = result.d_.val()
    grad_val(i) = result.d_.val();

    // ll only needs to be captured once (same value for every sweep)
#pragma omp critical
    {
      if (!ll_captured) { ll_val = result.val_.val(); ll_captured = true; }
    }

    // Reverse sweep on the inner var tape to get row i of the Hessian
    stan::math::grad(result.d_.vi_);

    for (int j = 0; j < p; ++j)
      H(i, j) = theta_fvd(j).val_.adj();
    // nested destructor recovers this thread's nested stack frame
  }
  } // end omp parallel (thread_stack destructor cleans up TLS)

  // Symmetrise (should already be symmetric to machine precision)
  H = 0.5 * (H + H.transpose());

  arma::vec g(grad_val.data(), p, true);
  arma::mat Hmat(H.data(), p, p, true);

  return Rcpp::List::create(
    Rcpp::Named("mean")     = ll_val,
    Rcpp::Named("gradient") = g,
    Rcpp::Named("Hessian")  = Hmat
  );
}
