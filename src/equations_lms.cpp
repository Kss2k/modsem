#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>
#include <vector>
#include <string>

#include "lms.h"
#include "pml.h"
#include "utils.h"
#include "mvnorm.h"

// [[Rcpp::depends( RcppArmadillo)]]


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
  const arma::mat Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

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
  // const arma::mat e = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else       zVec = arma::zeros<arma::vec>(numXis);

  const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);
  const arma::mat Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));

  const arma::mat Sxx = lX * A * Oi * A.t() * lX.t();
  const arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
  const arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
  const arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() +
    lY * (Binv * Psi * Binv.t()) * lY.t();

  return arma::join_cols(
    arma::join_rows(Sxx, Sxy),
    arma::join_rows(Sxy.t(), Syy)
  ) + d;
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



// Helper (assumed available in your TU)
extern arma::vec make_zvec(unsigned k, unsigned numXis, const arma::vec& z);
extern arma::mat make_Oi(unsigned k, unsigned numXis);


struct LMSModel {
  arma::mat A, Oxx, Oex, Ie, lY, lX, tY, tX, Gx, Ge,
    a, beta0, Psi, d, e, thr;

  arma::uvec isOrderedEnum; // if 0 variable is continuous, else idx in thesholds

  unsigned  k       = 0;
  unsigned  numXis  = 0;

  double pml;

  explicit LMSModel(const Rcpp::List& modFilled) {
    Rcpp::List matrices = modFilled["matrices"];
    Rcpp::List info     = modFilled["info"];
    Rcpp::List quad     = modFilled["quad"];

    k       = Rcpp::as<unsigned>(quad["k"]);
    numXis  = Rcpp::as<unsigned>(info["numXis"]);
    
    pml = Rcpp::as<std::string>(info["estimator"]) == "pml";

    A       = Rcpp::as<arma::mat>(matrices["A"]);
    Oxx     = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
    Oex     = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
    Ie      = Rcpp::as<arma::mat>(matrices["Ieta"]);
    lY      = Rcpp::as<arma::mat>(matrices["lambdaY"]);
    lX      = Rcpp::as<arma::mat>(matrices["lambdaX"]);
    tY      = Rcpp::as<arma::mat>(matrices["tauY"]);
    tX      = Rcpp::as<arma::mat>(matrices["tauX"]);
    Gx      = Rcpp::as<arma::mat>(matrices["gammaXi"]);
    Ge      = Rcpp::as<arma::mat>(matrices["gammaEta"]);
    a       = Rcpp::as<arma::mat>(matrices["alpha"]);
    beta0   = Rcpp::as<arma::mat>(matrices["beta0"]);
    Psi     = Rcpp::as<arma::mat>(matrices["psi"]);
    d       = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
    e       = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);
    thr     = Rcpp::as<arma::mat>(matrices["thresholds"]);

    isOrderedEnum = Rcpp::as<arma::uvec>(info["isOrderedEnum"]);

  }

  arma::vec mu(const arma::vec& z) const {
    const arma::vec zVec = make_zvec(k, numXis, z);
    const arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);
    const arma::mat Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

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
    const arma::mat Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    const arma::mat Oi = make_Oi(k, numXis);
    const arma::mat Sxx = lX * A * Oi * A.t() * lX.t();
    const arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
    const arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
    const arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() +
      lY * (Binv * Psi * Binv.t()) * lY.t();

    return arma::join_cols(
        arma::join_rows(Sxx, Sxy),
        arma::join_rows(Sxy.t(), Syy)
        ) + d;
  }

  LMSModel thread_clone() const {
    LMSModel c = *this;    // shallow copy baseline

    // Deep-copy matrices that can be modified downstream
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
    c.thr   = arma::mat(thr);

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
    case 14: return  M.thr  (r,c);
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
    const arma::uvec& n,
    const arma::uvec& d,
    const int npatterns) {
  const std::size_t J = V.n_rows;
  double ll = 0.0;

  for (std::size_t j = 0; j < J; j++) {
    if (arma::sum(TGamma[j]) <= DBL_MIN) continue;

    const arma::vec z  = V.row(j).t();   // view
    const arma::vec mu = M.mu(z);
    const arma::mat Sig = M.Sigma(z);

    for (int i = 0; i < npatterns; i++) {
      const arma::vec& nu = MeanPatterns[j][i];
      const arma::mat& S  = CovPatterns[j][i];
      const double tg     = TGamma[j][i];
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
                            const int npatterns = 1,
                            const int ncores = 1) {
  const LMSModel model = LMSModel(modelR);

  const arma::mat V       = Rcpp::as<arma::mat>(P["V"]);
  const arma::mat Pmat    = Rcpp::as<arma::mat>(P["P"]);
  const auto TGamma       = as_vec_of_vec(P["tgamma"]);
  const auto Mean         = as_vec_of_vec_of_vec(P["mean"]);
  const auto Cov          = as_vec_of_vec_of_mat(P["cov"]);
  const auto colidx       = as_vec_of_uvec(colidxR);

  return completeLogLikFromModel(
      model, V, TGamma,
      Mean, Cov, colidx, n,
      d, npatterns
  );
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
  const arma::mat Pmat    = Rcpp::as<arma::mat>(P["P"]);
  const auto      TGamma  = as_vec_of_vec(P["tgamma"]);
  const auto      Mean    = as_vec_of_vec_of_vec(P["mean"]);
  const auto      Cov     = as_vec_of_vec_of_mat(P["cov"]);
  const auto      colidx  = as_vec_of_uvec(colidxR);
  const Rcpp::List info   = modelR["info"];

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(
        mod, V, TGamma,
        Mean, Cov, colidx, n,
        d, npatterns
    );
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
                                const Rcpp::List& dataR,
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
  const arma::mat  Pmat    = Rcpp::as<arma::mat>(P["P"]);
  const auto       TGamma  = as_vec_of_vec(P["tgamma"]);
  const auto       Mean    = as_vec_of_vec_of_vec(P["mean"]);
  const auto       Cov     = as_vec_of_vec_of_mat(P["cov"]);
  const auto       colidx  = as_vec_of_uvec(colidxR);
  const auto       data    = as_vec_of_mat(dataR);

  const Rcpp::List info   = modelR["info"];

  auto comp_ll = [&](LMSModel& mod) -> double {
    return completeLogLikFromModel(
        mod, V, TGamma,
        Mean, Cov, colidx, n,
        d, npatterns
    );
  };

  return fdHessCpp(M, comp_ll, block, row, col, symmetric,
      relStep, minAbs, ncores); // multi-threaded
}


// stable elementwise log-sum-exp for two equally-sized matrices
inline arma::mat elem_logsumexp(const arma::mat& A, const arma::mat& B) {
  arma::mat M = arma::max(A, B);
  return M + arma::log(arma::exp(A - M) + arma::exp(B - M));
}


// row-wise log-sum-exp for a (L x Q) matrix
inline arma::vec row_logsumexp(const arma::mat& M) {
  arma::vec m = arma::max(M, 1);                     // L x 1
  arma::mat tmp = M.each_col() - m;
  return m + arma::log(arma::sum(arma::exp(tmp), 1));
}


// safe log of probabilities (clamps to keep finite)
inline arma::vec safe_log_vec(const arma::vec& v, double eps = 1e-300) {
  arma::vec out(v.n_elem, arma::fill::value(std::log(eps)));
  arma::uvec ok = arma::find(v > 0.0);
  if (!ok.is_empty()) out.elem(ok) = arma::log(v.elem(ok));
  return out;
}


// ---------- Pair-count container (per pattern, per pair) ----------
struct PairCount {
  arma::uword r_local, s_local;  // local indices (0-based) in this pattern
  arma::uword K_r, K_s;
  arma::umat  C;                 // K_r x K_s counts
  arma::uvec  rgrid, sgrid;      // flattened (column-major) 1..K indices, length K_r*K_s
};

struct PatternPCS {
  std::vector<PairCount> pcs;  // precomputed pair counts and grids
  arma::uvec isOrd_sub;        // per-pattern (copied from model$isOrderedEnum[idx])
  arma::uvec idx;              // columns indices into full (X,Y)
  arma::uword p_j{0}, n_j{0};
};

struct PCSBundle {
  std::vector<PatternPCS> patterns;  // size = npatterns
};


// ---- build counts for one ordinal pattern ----
static inline std::vector<PairCount>
precompute_pair_counts(const arma::mat& Xj,
                       const arma::uvec& isOrd_sub,
                       const arma::mat& thresholds) {
  const arma::uword p_j = Xj.n_cols, n_j = Xj.n_rows;
  std::vector<PairCount> out;
  if (p_j < 2u || n_j == 0u) return out;

  // If all rows have same width (common in practice):
  const arma::uword K_common = (thresholds.n_cols > 0 ? thresholds.n_cols - 1u : 0u);

  auto get_K = [&](arma::uword ord_row_1based) -> arma::uword {
    (void)ord_row_1based;
    return K_common;
    // If thresholds rows are ragged, replace with:
    // const arma::rowvec row = thresholds.row(ord_row_1based - 1u);
    // arma::uword K = row.n_elem - 1u;
    // return K;
  };

  for (arma::uword r = 0; r < p_j; ++r) {
    for (arma::uword s = r + 1u; s < p_j; ++s) {
      const arma::uword Kr = get_K(isOrd_sub[r]);
      const arma::uword Ks = get_K(isOrd_sub[s]);
      arma::umat C(Kr, Ks, arma::fill::zeros);

      arma::uvec rcat = arma::conv_to<arma::uvec>::from(Xj.col(r));
      arma::uvec scat = arma::conv_to<arma::uvec>::from(Xj.col(s));
      for (arma::uword t = 0; t < n_j; ++t) {
        arma::uword rr = rcat[t], ss = scat[t];
        if (rr >= 1u && rr <= Kr && ss >= 1u && ss <= Ks) C(rr-1u, ss-1u)++;
      }

      arma::uvec rgrid(Kr*Ks), sgrid(Kr*Ks);
      arma::uword k = 0;
      for (arma::uword ss = 1; ss <= Ks; ++ss)
        for (arma::uword rr = 1; rr <= Kr; ++rr, ++k) {
          rgrid[k] = rr; sgrid[k] = ss;
        }

      PairCount pc;
      pc.r_local = r; pc.s_local = s;
      pc.K_r = Kr; pc.K_s = Ks;
      pc.C = std::move(C);
      pc.rgrid = std::move(rgrid);
      pc.sgrid = std::move(sgrid);
      out.emplace_back(std::move(pc));
    }
  }
  return out;
}

// ---- build all patterns into a bundle and export as XPtr ----
// [[Rcpp::export]]
SEXP buildPCS_Xptr(const Rcpp::List& dataR,          // data$data.split
                   const Rcpp::List& colidxR,        // data$colidx
                   const arma::uvec& isOrderedEnum,  // model$info$isOrderedEnum
                   const arma::mat& thresholds)      // model$matrices$thresholds
{
  const auto data   = as_vec_of_mat(dataR);
  const auto colidx = as_vec_of_uvec(colidxR);

  if (data.size() != colidx.size())
    Rcpp::stop("buildPCS_Xptr: data and colidx lengths differ.");

  auto* bundle = new PCSBundle();
  bundle->patterns.resize(data.size());

  for (std::size_t j = 0; j < data.size(); ++j) {
    PatternPCS pj;
    pj.idx       = colidx[j];
    pj.isOrd_sub = isOrderedEnum.elem(pj.idx);
    pj.p_j = data[j].n_cols;
    pj.n_j = data[j].n_rows;

    const bool all_ordinal = arma::all(pj.isOrd_sub != 0u);
    if (all_ordinal && pj.p_j >= 2u && pj.n_j > 0u) {
      pj.pcs = precompute_pair_counts(data[j], pj.isOrd_sub, thresholds);
    } // else leave pcs empty -> fall back to general path

    bundle->patterns[j] = std::move(pj);
  }

  Rcpp::XPtr<PCSBundle> xptr(bundle, true);
  return xptr;
}


// Compute pairwise composite *log*-likelihood for one component from counts
// - mu_sub, Sig_sub: restricted to this pattern's variables
// - isOrd_sub: length p_j (>0 for all)
// - thresholds: full threshold matrix
// - pcs: precomputed counts/grids for this pattern
static inline double
pairwise_LL_from_counts_single(
    const arma::vec& mu_sub,
    const arma::mat& Sig_sub,
    const arma::uvec& isOrd_sub,
    const arma::mat& thresholds,
    const std::vector<PairCount>& pcs)
{
  double ll = 0.0;
  for (const auto& pc : pcs) {
    const arma::uword r = pc.r_local, s = pc.s_local;

    // Scalars as length-1 vectors for broadcasting in foo_vec_arma
    arma::vec mi(1), mj(1), Sii(1), Sjj(1), Sij(1);
    mi[0] = mu_sub[r];  mj[0] = mu_sub[s];
    Sii[0] = Sig_sub(r,r);  Sjj[0] = Sig_sub(s,s);  Sij[0] = Sig_sub(r,s);

    // thresholds rows
    const arma::vec tau_r = thresholds.row(isOrd_sub[r]-1u).t();
    const arma::vec tau_s = thresholds.row(isOrd_sub[s]-1u).t();

    // Probabilities for all cells (vector length = K_r*K_s), via one batched call
    arma::vec p_cells = foo_vec_arma(pc.rgrid, pc.sgrid, mi, mj, Sii, Sjj, Sij, tau_r, tau_s);

    // Accumulate: sum_{a,b} C[a,b] * log p[a,b]
    // We flattened C in column-major order to match rgrid/sgrid.
    arma::vec counts_vec = arma::conv_to<arma::vec>::from(arma::vectorise(pc.C)); // same order
    // stable/safe log to avoid -Inf * 0; your safe_log_fast works too
    arma::vec logp = arma::log(arma::clamp(p_cells, 1e-300, 1.0));
    ll += arma::dot(counts_vec, logp);
  }
  return ll;
}


inline arma::vec safe_log_local(const arma::vec& v, double eps = 1e-300) {
  arma::vec out(v.n_elem, arma::fill::value(std::log(eps)));
  arma::uvec ok = arma::find(v > 0.0);
  if (!ok.is_empty()) out.elem(ok) = arma::log(v.elem(ok));
  return out;
}


// Returns (n_j x npairs_j) matrix of *log* pairwise contributions for one component.
// - Continuous–Continuous (CC): log bivariate normal density
// - Ordinal–Continuous (OC/CO): log [ φ(x_cont) * Pr(interval | x_cont) ]
// - Ordinal–Ordinal (OO): log rectangle probability
inline arma::mat log_bvns_for_component(
    const arma::mat& Dj,            // n_j x p_j data for pattern j
    const arma::vec& mu_sub,        // length p_j
    const arma::mat& Sig_sub,       // p_j x p_j
    const arma::uvec& isOrd_sub,    // length p_j (0=cont, >0=row index in thresholds)
    const arma::mat& thresholds,    // thresholds matrix (rows indexed by isOrderedEnum-1)
    const int /*ncores*/ = 1
) {
  const arma::uword n_j = Dj.n_rows;
  const arma::uword p_j = Dj.n_cols;

  if (p_j < 2u || n_j == 0u)
    return arma::mat(n_j, 0u, arma::fill::zeros);

  const arma::uword npairs = p_j * (p_j - 1u) / 2u;
  arma::mat out(n_j, npairs);

  // Split columns once (so we can reuse categorized vs continuous storage)
  std::vector<arma::vec>  ccols(p_j);
  std::vector<arma::uvec> ocols(p_j);
  for (arma::uword j = 0; j < p_j; ++j) {
    if (isOrd_sub[j]) {
      ocols[j] = arma::conv_to<arma::uvec>::from(Dj.col(j));
    } else {
      ccols[j] = Dj.col(j);
    }
  }

  arma::uword col = 0u;
  for (arma::uword r = 0; r < p_j; ++r) {
    for (arma::uword s = r + 1u; s < p_j; ++s, ++col) {
      const bool ord_r = (isOrd_sub[r] != 0u);
      const bool ord_s = (isOrd_sub[s] != 0u);

      // Pull scalars/vectors in the shape expected by your kernels
      arma::vec mi(1), mj(1), Sii(1), Sjj(1), Sij(1);
      mi[0] = mu_sub[r]; mj[0] = mu_sub[s];
      Sii[0] = Sig_sub(r, r); Sjj[0] = Sig_sub(s, s); Sij[0] = Sig_sub(r, s);

      arma::vec prob_or_dens(n_j, arma::fill::zeros);

      if (!ord_r && !ord_s) {
        // CC
        const arma::vec& xr = ccols[r];
        const arma::vec& xs = ccols[s];
        prob_or_dens = fcc_vec_arma(xr, xs, mi, mj, Sii, Sjj, Sij); // density
      }
      else if (ord_r && !ord_s) {
        // OC: (ordinal r, continuous s) -> foc(x_cont = xs, r_ord = r, ...; tau_r)
        const arma::uvec& rr = ocols[r];
        const arma::vec&  xs = ccols[s];
        const arma::vec tau_r = thresholds.row(isOrd_sub[r] - 1u).t();
        // foc: (xj, r, mj, mk, Sjj, Skk, Sjk, tau_k)
        // map: xj=xs, r=rr, mj=mj (mean of xj/cont s), mk=mi (mean of ord r),
        //      Sjj=Sjj(s,s), Skk=Sii(r,r), Sjk=Sij(r,s)
        prob_or_dens = foc_vec_arma(xs, rr, mj, mi, Sjj, Sii, Sij, tau_r);
      }
      else if (!ord_r && ord_s) {
        // CO: symmetric -> foc(x_cont = xr, r_ord = s; tau_s)
        const arma::vec&  xr = ccols[r];
        const arma::uvec& ss = ocols[s];
        const arma::vec tau_s = thresholds.row(isOrd_sub[s] - 1u).t();
        // map: xj=xr, r=ss, mj=mi, mk=mj, Sjj=Sii, Skk=Sjj, Sjk=Sij
        prob_or_dens = foc_vec_arma(xr, ss, mi, mj, Sii, Sjj, Sij, tau_s);
      }
      else {
        // OO
        const arma::uvec& rr = ocols[r];
        const arma::uvec& ss = ocols[s];
        const arma::vec tau_r = thresholds.row(isOrd_sub[r] - 1u).t();
        const arma::vec tau_s = thresholds.row(isOrd_sub[s] - 1u).t();
        prob_or_dens = foo_vec_arma(rr, ss, mi, mj, Sii, Sjj, Sij, tau_r, tau_s); // probability
      }

      // Store *log* contribution, safely (0 -> log(eps))
      out.col(col) = safe_log_local(prob_or_dens);
    }
  }

  return out;  // n_j x npairs (LOG scale)
}



// Mixture sits *inside* each pair: for each pair do log-sum-exp over components,
// then sum across pairs and observations.
static inline double observedLogLikFromModelPML(
    const LMSModel&  M,
    const arma::mat& V,
    const arma::vec& w,
    const std::vector<arma::mat>& data,
    const std::vector<arma::uvec>& colidx,
    const arma::uvec& n,
    const PCSBundle* pcsBundle,       // may be nullptr
    const int npatterns = 1,
    const int ncores    = 1)
{
  const std::size_t Q = V.n_rows;
  double ll = 0.0;

  for (int j = 0; j < npatterns; ++j) {
    const arma::mat& Dj   = data[j];
    const arma::uvec& idx = colidx[j];
    const arma::uvec  isOrd_sub = M.isOrderedEnum.elem(idx);
    const arma::uword n_j = Dj.n_rows, p_j = Dj.n_cols;
    if (p_j < 2u || n_j == 0u) continue;

    // Try fast all-ordinal path via precomputed counts
    const bool all_ordinal = arma::all(isOrd_sub != 0u);
    const bool have_pcs =
      (pcsBundle != nullptr) &&
      (j < static_cast<int>(pcsBundle->patterns.size())) &&
      (!pcsBundle->patterns[j].pcs.empty());

    if (all_ordinal && have_pcs) {
      const auto& pcs = pcsBundle->patterns[j].pcs;

      // Precompute per-q restricted moments
      std::vector<arma::vec> mu_sub_Q(Q);
      std::vector<arma::mat> Sig_sub_Q(Q);
      for (std::size_t q = 0; q < Q; ++q) {
        if (w[q] <= DBL_MIN) continue;
        const arma::vec z   = V.row(q).t();
        const arma::vec mu  = M.mu(z);
        const arma::mat Sig = M.Sigma(z);
        mu_sub_Q[q]  = mu.elem(idx);
        Sig_sub_Q[q] = Sig.submat(idx, idx);
      }

      // For each pair, do mixture inside pair at the cell level
      for (const auto& pc : pcs) {
        const arma::uword r = pc.r_local, s = pc.s_local;
        const arma::uword L = pc.rgrid.n_elem;

        arma::mat LW(L, Q, arma::fill::value(-std::numeric_limits<double>::infinity()));
        for (std::size_t q = 0; q < Q; ++q) {
          if (w[q] <= DBL_MIN) continue;

          const arma::vec& mu_sub  = mu_sub_Q[q];
          const arma::mat& Sig_sub = Sig_sub_Q[q];

          arma::vec mi(1), mj(1), Sii(1), Sjj(1), Sij(1);
          mi[0] = mu_sub[r];  mj[0] = mu_sub[s];
          Sii[0] = Sig_sub(r,r);  Sjj[0] = Sig_sub(s,s);  Sij[0] = Sig_sub(r,s);

          const arma::vec tau_r = M.thr.row(isOrd_sub[r]-1u).t();
          const arma::vec tau_s = M.thr.row(isOrd_sub[s]-1u).t();

          arma::vec pq = foo_vec_arma(pc.rgrid, pc.sgrid, mi, mj, Sii, Sjj, Sij, tau_r, tau_s);
          LW.col(q) = std::log(std::max(w[q], DBL_MIN)) + safe_log_vec(pq);
        }

        arma::vec logmix = row_logsumexp(LW);
        arma::vec counts_vec = arma::conv_to<arma::vec>::from(arma::vectorise(pc.C));
        ll += arma::dot(counts_vec, logmix);
      }
      continue;
    }

    // General path (mixed/continuous) using per-row pairwise logs
    bool first = true;
    arma::mat logmix_j;
    for (std::size_t q = 0; q < Q; ++q) {
      if (w[q] <= DBL_MIN) continue;
      const arma::vec z   = V.row(q).t();
      const arma::vec mu  = M.mu(z);
      const arma::mat Sig = M.Sigma(z);
      arma::vec  mu_sub  = mu.elem(idx);
      arma::mat  Sig_sub = Sig.submat(idx, idx);

      arma::mat Lq = log_bvns_for_component(Dj, mu_sub, Sig_sub, isOrd_sub, M.thr, ncores);
      Lq += std::log(std::max(w[q], DBL_MIN));

      if (first) { logmix_j = Lq; first = false; }
      else       { logmix_j = elem_logsumexp(logmix_j, Lq); }
    }
    if (!first) {
      ll += arma::sum(arma::sum(logmix_j, 1));
    }
  }
  return ll;
}


// [[Rcpp::export]]
double observedLogLikLmsPMLCpp(const Rcpp::List& modelR,
    const Rcpp::List& dataR, const Rcpp::List& colidxR, const Rcpp::List& P,
    SEXP pcs_xptr,                      // external pointer from buildPCS_Xptr()
    const arma::uvec& n, const int npatterns = 1L, const int ncores = 1L) {
  const LMSModel M(modelR);
  const arma::mat V = Rcpp::as<arma::mat>(P["V"]);
  const arma::vec w = Rcpp::as<arma::vec>(P["w"]);
  const auto colidx = as_vec_of_uvec(colidxR);
  const auto data   = as_vec_of_mat(dataR);

  const PCSBundle* pcsBundle =
      (pcs_xptr == R_NilValue) ? nullptr : Rcpp::XPtr<PCSBundle>(pcs_xptr).get();

  return observedLogLikFromModelPML(M, V, w, data, colidx, n, pcsBundle, npatterns, ncores);
}


// [[Rcpp::export]]
arma::vec gradObsLogLikLmsPMLCpp(const Rcpp::List& modelR,
                                     const Rcpp::List& dataR,
                                     const Rcpp::List& colidxR,
                                     const Rcpp::List& P,
                                     SEXP pcs_xptr,
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
  const PCSBundle* pcsBundle =
      (pcs_xptr == R_NilValue) ? nullptr : Rcpp::XPtr<PCSBundle>(pcs_xptr).get();

  auto obs_ll = [&](LMSModel& mod) -> double {
    return observedLogLikFromModelPML(mod, V, w, data, colidx, n, pcsBundle, npatterns, 1L);
  };
  return gradientFD(M, obs_ll, block, row, col, symmetric, eps, ncores);
}


// [[Rcpp::export]]
Rcpp::List hessObsLogLikLmsPMLCpp(const Rcpp::List& modelR,
                                      const Rcpp::List& dataR,
                                      const Rcpp::List& colidxR,
                                      const Rcpp::List& P,
                                      SEXP pcs_xptr,
                                      const arma::uvec& block,
                                      const arma::uvec& row,
                                      const arma::uvec& col,
                                      const arma::uvec& symmetric,
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
  const PCSBundle* pcsBundle =
      (pcs_xptr == R_NilValue) ? nullptr : Rcpp::XPtr<PCSBundle>(pcs_xptr).get();

  auto obs_ll = [&](LMSModel& mod) -> double {
    return observedLogLikFromModelPML(mod, V, w, data, colidx, n, pcsBundle, npatterns, 1L);
  };
  return fdHessCpp(M, obs_ll, block, row, col, symmetric, relStep, minAbs, ncores);
}
