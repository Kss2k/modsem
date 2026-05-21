#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>

#include "lms.h"
#include "utils.h"
#include "mvnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]


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
  arma::mat A, Oxx, Oex, Ie, lY, lX, W, T, tY, tX, Gx, Ge,
    a, beta0, Psi, d, e;
  unsigned  k        = 0;
  unsigned  numXis   = 0;
  bool hasComposites = false;

  // Precomputed z-independent derived quantities; call update_cache() after any
  // change to A or Gx (i.e. after lms_param / setParams perturbations).
  arma::mat Oi, GxA, AOi, AOiAt;

  void update_cache() {
    Oi    = make_Oi(k, numXis);
    GxA   = Gx * A;
    AOi   = A * Oi;
    AOiAt = AOi * A.t();
  }

  explicit LMSModel(const Rcpp::List& modFilled) {

    Rcpp::List matrices = modFilled["matrices"];
    Rcpp::List info     = modFilled["info"];
    Rcpp::List quad     = modFilled["quad"];

    k       = Rcpp::as<unsigned>(quad["k"]);
    numXis  = Rcpp::as<unsigned>(info["numXis"]);
    hasComposites = Rcpp::as<bool>(info["hasComposites"]);

    // one-liners, no loops
    A      = Rcpp::as<arma::mat>(matrices["A"]);
    Oxx    = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
    Oex    = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
    Ie     = Rcpp::as<arma::mat>(matrices["Ieta"]);
    lY     = Rcpp::as<arma::mat>(matrices["lambdaY"]);
    lX     = Rcpp::as<arma::mat>(matrices["lambdaX"]);
    tY     = Rcpp::as<arma::mat>(matrices["tauY"]);
    tX     = Rcpp::as<arma::mat>(matrices["tauX"]);
    W      = Rcpp::as<arma::mat>(matrices["W"]);
    T      = Rcpp::as<arma::mat>(matrices["T"]);
    Gx     = Rcpp::as<arma::mat>(matrices["gammaXi"]);
    Ge     = Rcpp::as<arma::mat>(matrices["gammaEta"]);
    a      = Rcpp::as<arma::mat>(matrices["alpha"]);
    beta0  = Rcpp::as<arma::mat>(matrices["beta0"]);
    Psi    = Rcpp::as<arma::mat>(matrices["psi"]);
    d      = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
    e      = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);

    update_cache();
  }

  // Combined mu+Sigma: avoids recomputing zVec/kronZ/Binv/lXc for the same z.
  // Use this in any hot path that needs both.
  std::pair<arma::vec, arma::mat> muSigma(const arma::vec& z) const {
    const arma::vec zVec  = make_zvec(k, numXis, z);
    const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);
    const arma::mat Binv  = arma::inv(Ie - Ge - kronZ.t() * Oex);

    const arma::vec muXi  = beta0 + A * zVec;
    const arma::vec muEta = Binv * (a +
        Gx * muXi +
        kronZ.t() * Oxx * muXi);

    const arma::mat Eta      = Binv * (GxA + kronZ.t() * Oxx * A);
    const arma::mat varXi    = AOiAt;
    const arma::mat varEta   = Eta * Oi * Eta.t() + Binv * Psi * Binv.t();
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
    c.W     = arma::mat(W);
    c.T     = arma::mat(T);
    c.Gx    = arma::mat(Gx);
    c.Ge    = arma::mat(Ge);
    c.a     = arma::mat(a);
    c.beta0 = arma::mat(beta0);
    c.Psi   = arma::mat(Psi);
    c.d     = arma::mat(d);
    c.e     = arma::mat(e);
    c.Oi    = arma::mat(Oi);
    c.GxA   = arma::mat(GxA);
    c.AOi   = arma::mat(AOi);
    c.AOiAt = arma::mat(AOiAt);

    return c;
  }
};


// [[Rcpp::export]]
Rcpp::List muSigmaLmsCpp(Rcpp::List model, arma::vec z) {
  const std::pair<arma::vec, arma::mat> ms = LMSModel(model).muSigma(z);
  return Rcpp::List::create(
    Rcpp::Named("mu")    = ms.first,
    Rcpp::Named("sigma") = ms.second
  );
}


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
    case 14: return  M.W    (r,c);
    case 15: return  M.T    (r,c);
    default: Rcpp::stop("unknown block id");
  }
}

struct LMSAdjoints {
  arma::mat A, Oxx, Oex, lX, tX, Gx, Ge, a, beta0, Psi, d;

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
    d(M.d.n_rows, M.d.n_cols, arma::fill::zeros) {}

  void add(const LMSAdjoints& other) {
    A     += other.A;
    Oxx   += other.Oxx;
    Oex   += other.Oex;
    lX    += other.lX;
    tX    += other.tX;
    Gx    += other.Gx;
    Ge    += other.Ge;
    a     += other.a;
    beta0 += other.beta0;
    Psi   += other.Psi;
    d     += other.d;
  }
};


inline void add_kron_Ie_u_adjoint(arma::vec& u_bar,
                                  const arma::mat& K_bar,
                                  const unsigned numXis,
                                  const unsigned numEtas) {
  for (unsigned e = 0; e < numEtas; ++e) {
    u_bar += K_bar.submat(e * numXis, e, (e + 1) * numXis - 1, e);
  }
}


inline const arma::mat& lms_adjoint_block(const LMSAdjoints& adj,
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
    default: {
      static const arma::mat empty;
      return empty;
    }
  }
}


inline void accumulateMuSigmaAdjoints(const LMSModel& M,
                                      const arma::vec& z,
                                      const arma::vec& mu_bar,
                                      const arma::mat& Sigma_bar,
                                      LMSAdjoints& adj) {
  const arma::uword numEta = M.Ie.n_rows;
  const arma::vec zVec = make_zvec(M.k, M.numXis, z);
  const arma::vec u = M.beta0 + M.A * zVec;
  const arma::mat K = arma::kron(M.Ie, u);
  const arma::mat B = M.Ie - M.Ge - K.t() * M.Oex;
  const arma::mat Binv = arma::inv(B);
  const arma::vec r = M.a + M.Gx * u + K.t() * M.Oxx * u;
  const arma::vec v = Binv * r;
  const arma::mat C = M.Gx * M.A + K.t() * M.Oxx * M.A;
  const arma::mat Eta = Binv * C;
  const arma::mat varXi = M.AOiAt;
  const arma::mat varEta = Eta * M.Oi * Eta.t() + Binv * M.Psi * Binv.t();
  const arma::mat covXiEta = M.AOi * Eta.t();
  const arma::mat V = arma::join_cols(
    arma::join_rows(varXi,        covXiEta),
    arma::join_rows(covXiEta.t(), varEta)
  );
  const arma::vec xieta = arma::join_cols(u, v);

  if (adj.tX.n_elem == mu_bar.n_elem) adj.tX += mu_bar;
  adj.lX += mu_bar * xieta.t();

  arma::vec xieta_bar = M.lX.t() * mu_bar;
  arma::mat V_bar = M.lX.t() * Sigma_bar * M.lX;
  adj.lX += (Sigma_bar + Sigma_bar.t()) * M.lX * V;
  adj.d += Sigma_bar;

  arma::vec u_bar = xieta_bar.subvec(0, M.numXis - 1);
  arma::vec v_bar = xieta_bar.subvec(M.numXis, M.numXis + numEta - 1);

  arma::mat varXi_bar = V_bar.submat(0, 0, M.numXis - 1, M.numXis - 1);
  arma::mat cov_bar = V_bar.submat(0, M.numXis, M.numXis - 1, M.numXis + numEta - 1);
  arma::mat varEta_bar = V_bar.submat(M.numXis, M.numXis,
                                      M.numXis + numEta - 1,
                                      M.numXis + numEta - 1);
  cov_bar += V_bar.submat(M.numXis, 0,
                          M.numXis + numEta - 1,
                          M.numXis - 1).t();

  adj.A += (varXi_bar + varXi_bar.t()) * M.AOi;

  arma::mat Eta_bar = cov_bar.t() * M.AOi;
  adj.A += cov_bar * Eta * M.Oi;

  Eta_bar += (varEta_bar + varEta_bar.t()) * Eta * M.Oi;
  arma::mat Binv_bar = varEta_bar * Binv * M.Psi.t() +
                       varEta_bar.t() * Binv * M.Psi;
  adj.Psi += Binv.t() * varEta_bar * Binv;

  Binv_bar += Eta_bar * C.t();
  arma::mat C_bar = Binv.t() * Eta_bar;

  Binv_bar += v_bar * r.t();
  arma::vec r_bar = Binv.t() * v_bar;

  arma::mat B_bar = -Binv.t() * Binv_bar * Binv.t();
  adj.Ge -= B_bar;
  arma::mat K_bar = -M.Oex * B_bar.t();
  adj.Oex += -K * B_bar;

  adj.a += r_bar;
  adj.Gx += r_bar * u.t();
  u_bar += M.Gx.t() * r_bar;

  arma::vec q = M.Oxx * u;
  K_bar += q * r_bar.t();
  arma::vec q_bar = K * r_bar;
  adj.Oxx += q_bar * u.t();
  u_bar += M.Oxx.t() * q_bar;

  adj.Gx += C_bar * M.A.t();
  adj.A += M.Gx.t() * C_bar;

  arma::mat M_bar = K * C_bar;
  K_bar += (M.Oxx * M.A) * C_bar.t();
  adj.Oxx += M_bar * M.A.t();
  adj.A += M.Oxx.t() * M_bar;

  add_kron_Ie_u_adjoint(u_bar, K_bar, M.numXis, numEta);

  adj.beta0 += u_bar;
  adj.A += u_bar * zVec.t();
}


inline bool completeLogLikScore(const arma::vec& mu,
                                const arma::mat& sigma,
                                const arma::vec& nu,
                                const arma::mat& S,
                                const double tgamma,
                                arma::vec& mu_bar,
                                arma::mat& Sigma_bar) {
  if (!sigma.is_finite()) return false;

  arma::mat L;
  if (!arma::chol(L, sigma, "lower")) return false;

  const arma::vec diff = nu - mu;
  const arma::vec inv_diff = arma::solve(
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

  mu_bar = tgamma * inv_diff;
  Sigma_bar = 0.5 * (Sinv * (S + tgamma * diff * diff.t()) * Sinv -
                     tgamma * Sinv);

  return true;
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
    const bool set_threads = true) {

  LMSAdjoints adj(M);
  const std::size_t J = V.n_rows;
  const arma::uword pObs = M.lX.n_rows;
  const std::size_t npar = block.n_elem;
  bool failed = false;

#ifdef _OPENMP
  int old_threads = omp_get_max_threads();
  if (set_threads) {
    if (ncores <= 0)
      Rcpp::stop("ncores must be positive");
    omp_set_num_threads(ncores);
  }
  const bool in_parallel = omp_in_parallel();
  const int nthreads = in_parallel ? 1 : omp_get_max_threads();
#else
  const int nthreads = 1;
#endif
  std::vector<LMSAdjoints> thread_adj;
  thread_adj.reserve(nthreads);
  for (int t = 0; t < nthreads; ++t)
    thread_adj.emplace_back(M);

#pragma omp parallel for default(none) if(!in_parallel) \
  shared(M, V, TGamma, MeanPatterns, CovPatterns, colidx, J, pObs, npatterns, \
         thread_adj, in_parallel) reduction(||:failed) schedule(static)
  for (std::size_t j = 0; j < J; j++) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    const arma::vec z = V.row(j).t();
    const std::pair<arma::vec, arma::mat> ms = M.muSigma(z);
    const arma::vec& mu = ms.first;
    const arma::mat& Sig = ms.second;

    arma::vec mu_bar_full(pObs);
    arma::mat Sig_bar_full(pObs, pObs);
    mu_bar_full.zeros();
    Sig_bar_full.zeros();

    for (int i = 0; i < npatterns; i++) {
      const double tg = TGamma[j][i];
      if (tg <= DBL_MIN) continue;

      arma::vec mu_bar_i;
      arma::mat Sig_bar_i;
      const arma::uvec& idx = colidx[i];

      if (!completeLogLikScore(mu.elem(idx),
                               Sig.submat(idx, idx),
                               MeanPatterns[j][i],
                               CovPatterns[j][i],
                               tg,
                               mu_bar_i,
                               Sig_bar_i)) {
        failed = true;
        continue;
      }

      mu_bar_full.elem(idx) += mu_bar_i;
      Sig_bar_full.submat(idx, idx) += Sig_bar_i;
    }

    accumulateMuSigmaAdjoints(M, z, mu_bar_full, Sig_bar_full, thread_adj[tid]);
  }

  if (failed) {
    arma::vec grad(npar);
    grad.fill(arma::datum::nan);
#ifdef _OPENMP
    if (set_threads) omp_set_num_threads(old_threads);
#endif
    return grad;
  }

  for (int t = 0; t < nthreads; ++t)
    adj.add(thread_adj[t]);

  arma::vec grad(npar, arma::fill::zeros);
  for (std::size_t k = 0; k < npar; ++k) {
    const arma::mat& A = lms_adjoint_block(adj, block[k]);
    if (A.is_empty()) continue;

    grad[k] = A(row[k], col[k]);
    if (symmetric[k] && row[k] != col[k])
      grad[k] += A(col[k], row[k]);
  }

#ifdef _OPENMP
  if (set_threads) omp_set_num_threads(old_threads);
#endif

  return grad;
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
    Mc.update_cache();

    // Evaluate on the perturbed *local* model
    const double f1 = logLik(Mc);

    // Gradient component
    grad[k] = (f1 - f0) / eps;

    // No need to restore: Mc is thread-local and will be destroyed here.
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

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(mod, V, TGamma, Mean, Cov, colidx, n, d, npatterns);
  };

  if (!M.hasComposites) {
    return completeGradientReverseFromModel(M, V, TGamma, Mean, Cov, colidx,
                                            block, row, col, symmetric,
                                            npatterns, ncores);
  }

  return gradientFD(M, comp_ll, block, row, col, symmetric, eps, ncores);
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

    const arma::vec z = V.row(i).t();
    const std::pair<arma::vec, arma::mat> ms = M.muSigma(z);
    const arma::vec& mu  = ms.first;
    const arma::mat& Sig = ms.second;

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
    return observedLogLikFromModel(mod, V, w, samplingWeights, data, colidx, n, npatterns, 1L);
  };

  return gradientFD(M, obs_ll, block, row, col, symmetric, eps, ncores);
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
    Mc.update_cache();
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
    Mc.update_cache();
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
                               const double      relStep   = 1e-6,
                               const double      minAbs    = 0.0,
                               const int         ncores    = 1L) {
  LMSModel M(modelR);

  const arma::mat V               = Rcpp::as<arma::mat>(P["V"]);
  const arma::vec w               = Rcpp::as<arma::vec>(P["w"]);
  const arma::vec samplingWeights = Rcpp::as<arma::vec>(P["sampling.weights"]);

  const auto colidx = as_vec_of_uvec(colidxR);
  const auto data   = as_vec_of_mat(dataR);

  auto obs_ll = [&](LMSModel& mod) -> double {
    return observedLogLikFromModel(mod, V, w, samplingWeights, data, colidx, n, npatterns, 1L);
  };

  return fdHessCpp(M, obs_ll, block, row, col, symmetric, relStep, minAbs, ncores);
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
  const double min_scale = (minAbs > 0.0) ? minAbs : 1.0;
  const arma::vec incr =
    arma::max(arma::abs(base), arma::vec(p).fill(min_scale)) * relStep;

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
    LMSModel Mc = M.thread_clone();
    arma::vec pars = base;
    pars[j] += incr[j];
    setParams(Mc, block, row, col, symmetric, pars);
    Mc.update_cache();

    const arma::vec grad_j = completeGradientReverseFromModel(
      Mc, V, TGamma, Mean, Cov, colidx, block, row, col, symmetric,
      npatterns, 1L, false
    );

    Hess.col(j) = (grad_j - grad0) / incr[j];
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

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(mod, V, TGamma, Mean, Cov, colidx, n, d, npatterns);
  };

  if (!M.hasComposites) {
    return fdHessFromCompleteGradient(M, V, TGamma, Mean, Cov, colidx, n, d,
                                      block, row, col, symmetric, npatterns,
                                      relStep, minAbs, ncores);
  }

  return fdHessCpp(M, comp_ll, block, row, col, symmetric, relStep, minAbs, ncores);
}


// [[Rcpp::export]]
arma::mat densityMatrixLmsCpp(const Rcpp::List& modelR,
                               const arma::mat&  V,
                               const Rcpp::List& dataR,
                               const Rcpp::List& colidxR,
                               const arma::uvec& n,
                               const arma::vec&  samplingWeights,
                               const int npatterns = 1) {
  const LMSModel M(modelR);
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);

  const std::size_t Q  = V.n_rows;
  const int         N  = (int)arma::sum(n);
  const bool        hasSW = samplingWeights.n_elem == (std::size_t)N;

  arma::mat out(N, Q, arma::fill::zeros);

  for (std::size_t i = 0; i < Q; ++i) {
    const arma::vec z   = V.row(i).t();
    const auto      ms  = M.muSigma(z);
    const arma::vec& mu  = ms.first;
    const arma::mat& Sig = ms.second;

    int offset = 0;
    for (int j = 0; j < npatterns; ++j) {
      const int end = offset + (int)n[j] - 1;
      out(arma::span(offset, end), i) =
        dmvnfast(data[j], mu.elem(colidx[j]), Sig.submat(colidx[j], colidx[j]),
                 false, 1, false);
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
                                const int npatterns = 1) {
  const auto data = as_vec_of_mat(dataR);
  const std::size_t Q = P.n_cols;

  Rcpp::List wMeans(Q), wCovs(Q), tGamma(Q);

  for (std::size_t i = 0; i < Q; ++i) {
    const arma::vec p = P.col(i);

    Rcpp::List wMeans_i(npatterns), wCovs_i(npatterns);
    arma::vec  tGamma_i(npatterns);

    int offset = 0;
    for (int j = 0; j < npatterns; ++j) {
      const int        end = offset + (int)n[j] - 1;
      const arma::vec   pj = p.subvec(offset, end);
      const arma::mat&  Dj = data[j];

      const double tg  = arma::sum(pj);
      arma::vec    wm  = Dj.t() * pj / tg;
      arma::mat    X   = Dj.each_row() - wm.t();
      arma::mat    cov = X.t() * (X.each_col() % pj);

      wMeans_i[j] = wm;
      wCovs_i[j]  = cov;
      tGamma_i[j] = tg;

      offset = end + 1;
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
