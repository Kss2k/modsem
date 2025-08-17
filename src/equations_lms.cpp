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
  const arma::mat lY = Rcpp::as<arma::mat>(matrices["lambdaY"]);
  const arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::mat tY = Rcpp::as<arma::mat>(matrices["tauY"]);
  const arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  const arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  const arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else       zVec = arma::zeros<arma::vec>(numXis);

  const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);
  }

  const arma::vec muX = tX + lX * (beta0 + A * zVec);
  const arma::vec muY = tY +
    lY * (Binv * (a +
          Gx * (beta0 + A * zVec) +
          kronZ.t() * Oxx * (beta0 + A * zVec)));

  return arma::join_cols(muX, muY);
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
  const arma::mat lY = Rcpp::as<arma::mat>(matrices["lambdaY"]);
  const arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::mat tY = Rcpp::as<arma::mat>(matrices["tauY"]);
  const arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  const arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  const arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);
  const arma::mat Psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat d = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  const arma::mat e = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else       zVec = arma::zeros<arma::vec>(numXis);

  const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);
  }

  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));

  const arma::mat Sxx = lX * A * Oi * A.t() * lX.t() + d;
  const arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
  const arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
  const arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() +
    lY * (Binv * Psi * Binv.t()) * lY.t() + e;

  return arma::join_cols(arma::join_rows(Sxx, Sxy), arma::join_rows(Sxy.t(), Syy));
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

    arma::mat Binv;
    if (Ie.n_cols == 1) Binv = arma::mat(Ie);
    else                Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    const arma::vec muX = tX + lX * (beta0 + A * zVec);
    const arma::vec muY = tY +
      lY * (Binv * (a +
            Gx * (beta0 + A * zVec) +
            kronZ.t() * Oxx * (beta0 + A * zVec)));

    return arma::join_cols(muX, muY);
  }

  arma::mat Sigma(const arma::vec& z) const {
    const arma::vec zVec  = make_zvec(k, numXis, z);
    const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

    arma::mat Binv;
    if (Ie.n_cols == 1) Binv = arma::mat(Ie);
    else                Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    const arma::mat Oi = make_Oi(k, numXis);
    const arma::mat Sxx = lX * A * Oi * A.t() * lX.t() + d;
    const arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
    const arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
    const arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() +
      lY * (Binv * Psi * Binv.t()) * lY.t() + e;

    return arma::join_cols(
        arma::join_rows(Sxx, Sxy),
        arma::join_rows(Sxy.t(), Syy)
        );
  }

  LMSModel thread_clone() const {
    LMSModel c = *this;    // shallow for everything (fast)
                           // Deep-copy ONLY what set_params()/lms_param can modify:
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

      ll += totalDmvnWeightedCpp(
        mu.elem(colidx[i]),
        Sig.submat(colidx[i], colidx[i]),
        nu, S, tg, n[i], d[i]);
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

  const Rcpp::List info   = modelR["info"];

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(mod, V, TGamma, Mean, Cov,
                                   colidx, n, d, npatterns);
  };

  return gradientFD(M, comp_ll, block, row, col, symmetric, eps, ncores);
}




inline double observedLogLikFromModel(const LMSModel&  M,
                                      const arma::mat& V,
                                      const arma::vec& w,
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

  return arma::sum(arma::log(density));
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
  const auto colidx = as_vec_of_uvec(colidxR);
  const auto data   = as_vec_of_mat(dataR);

  auto obs_ll = [&](LMSModel& mod) -> double {
    return observedLogLikFromModel(mod, V, w, data, colidx, n, npatterns, 1L); // single-threaded
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
  const auto colidx = as_vec_of_uvec(colidxR);
  const auto data   = as_vec_of_mat(dataR);

  return observedLogLikFromModel(M, V, w, data, colidx, n, npatterns, ncores);
}


inline arma::vec get_params(const LMSModel& M,
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


inline void set_params(LMSModel&         M,
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

  // Baseline parameter vector and FD step sizes
  const arma::vec base = get_params(M, block, row, col);
  const arma::vec incr =
    arma::max(arma::abs(base), arma::vec(p).fill(minAbsPar)) * relStep;

  // Build Koschal displacement matrix (read-only afterwards)
  std::vector< arma::vec > disp;
  disp.reserve(1 + 2*p + (p*(p-1))/2); // rough lower bound
  disp.emplace_back(arma::zeros<arma::vec>(p));          // origin
  for (std::size_t i = 0; i < p; ++i) {                  // +e_i / –e_i
    arma::vec v = arma::zeros<arma::vec>(p);
    v[i] = 1;  disp.push_back(v);
    v[i] = -1; disp.push_back(v);
  }
  for (std::size_t i = 0; i < p - 1; ++i)                // +e_i + e_j (i<j)
    for (std::size_t j = i + 1; j < p; ++j) {
      arma::vec v = arma::zeros<arma::vec>(p);
      v[i] = v[j] = 1;
      disp.push_back(v);
    }
  const std::size_t m = disp.size();

  // Evaluate fun at every design point (parallel)
  arma::vec y(m);

#pragma omp parallel for default(none) \
  shared(M, disp, m, block, row, col, symmetric, base, incr, y) \
  firstprivate(fun) schedule(static)
  for (std::size_t k = 0; k < m; ++k) {
    // If you have LMSModel::shallow_clone(), prefer it:
    LMSModel Mc = M.thread_clone();

    // θ = base + disp[k] % incr
    set_params(Mc, block, row, col, symmetric, base + disp[k] % incr);

    // Evaluate on local model
    y[k] = fun(Mc);
    // No restore needed (Mc is thread-local)
  }

  // Restore θ₀ on the original model (serial, for callers expecting M unchanged)
  set_params(M, block, row, col, symmetric, base);

  // Build design matrix X (serial; cheap vs evals; BLAS may thread this)
  const std::size_t q = 1 + 2*p + (p*(p-1))/2;
  arma::mat X(m, q, arma::fill::ones);
  std::size_t col_id = 1;

  // linear terms
  for (std::size_t j = 0; j < p; ++j, ++col_id)
    for (std::size_t k = 0; k < m; ++k)
      X(k, col_id) = disp[k][j];

  // squares
  for (std::size_t j = 0; j < p; ++j, ++col_id)
    for (std::size_t k = 0; k < m; ++k)
      X(k, col_id) = std::pow(disp[k][j], 2);

  // cross terms
  for (std::size_t i = 0; i < p - 1; ++i)
    for (std::size_t j = i + 1; j < p; ++j, ++col_id)
      for (std::size_t k = 0; k < m; ++k)
        X(k, col_id) = disp[k][i] * disp[k][j];

  // “frac” scaling (nlme-compatible)
  arma::vec frac(q, arma::fill::ones);
  for (std::size_t j = 0; j < p; ++j)              frac[1 + j]     = incr[j];
  for (std::size_t j = 0; j < p; ++j)              frac[1 + p + j] = incr[j] * incr[j];
  col_id = 1 + 2*p;
  for (std::size_t i = 0; i < p - 1; ++i)
    for (std::size_t j = i + 1; j < p; ++j, ++col_id)
      frac[col_id] = incr[i] * incr[j];

  // Solve for polynomial coefficients
  arma::vec coef = arma::solve(X, y) / frac;

  // Gradient (first-order coefs)
  arma::vec grad = coef.subvec(1, p);

  // Hessian
  arma::mat Hess(p, p, arma::fill::zeros);
  for (std::size_t j = 0; j < p; ++j)               // diagonal: 2 * c_j
    Hess(j, j) = 2.0 * coef[1 + p + j];

  col_id = 1 + 2*p;                                 // off-diagonal: d_ij
  for (std::size_t i = 0; i < p - 1; ++i)
    for (std::size_t j = i + 1; j < p; ++j, ++col_id) {
      Hess(i, j) = coef[col_id];
      Hess(j, i) = coef[col_id];
    }

  return Rcpp::List::create(
      Rcpp::Named("mean")     = coef[0],
      Rcpp::Named("gradient") = grad,
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
                               const double      relStep = 1e-6,
                               const double      minAbs  = 0.0,
                               const int         ncores  = 1L) {
    LMSModel M(modelR);

    const arma::mat V = Rcpp::as<arma::mat>(P["V"]);
    const arma::vec w = Rcpp::as<arma::vec>(P["w"]);
    const auto colidx = as_vec_of_uvec(colidxR);
    const auto data   = as_vec_of_mat(dataR);

    auto obs_ll = [&](LMSModel& mod) -> double {
        return observedLogLikFromModel(mod, V, w, data, colidx,
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
