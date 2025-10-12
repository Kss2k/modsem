#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <float.h>
#include <cmath>

#include "lms.h"
#include "utils.h"
#include "mvnorm.h"
#include "lms_autodiff.h"

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

namespace {

using modsem::lms::Model;
using modsem::lms::CompleteData;
using modsem::lms::ObservedData;

inline Eigen::MatrixXd as_matrix(const Rcpp::RObject& obj) {
  return Rcpp::as<Eigen::MatrixXd>(obj);
}

inline Eigen::VectorXd as_vector(const Rcpp::RObject& obj) {
  Rcpp::NumericVector v(obj);
  Eigen::VectorXd res(v.size());
  for (R_xlen_t i = 0; i < v.size(); ++i) res[i] = v[i];
  return res;
}

inline Eigen::VectorXi to_eigen(const arma::uvec& x) {
  Eigen::VectorXi res(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) res[i] = static_cast<int>(x[i]);
  return res;
}

Model<double> build_model_eigen(const Rcpp::List& modelR) {
  Model<double> M;

  const Rcpp::List matrices = modelR["matrices"];
  const Rcpp::List info     = modelR["info"];
  const Rcpp::List quad     = modelR["quad"];

  M.k      = Rcpp::as<unsigned>(quad["k"]);
  M.numXis = Rcpp::as<unsigned>(info["numXis"]);

  M.A     = as_matrix(matrices["A"]);
  M.Oxx   = as_matrix(matrices["omegaXiXi"]);
  M.Oex   = as_matrix(matrices["omegaEtaXi"]);
  M.Ie    = as_matrix(matrices["Ieta"]);
  M.lY    = as_matrix(matrices["lambdaY"]);
  M.lX    = as_matrix(matrices["lambdaX"]);
  M.tY    = as_vector(matrices["tauY"]);
  M.tX    = as_vector(matrices["tauX"]);
  M.Gx    = as_matrix(matrices["gammaXi"]);
  M.Ge    = as_matrix(matrices["gammaEta"]);
  M.a     = as_vector(matrices["alpha"]);
  M.beta0 = as_vector(matrices["beta0"]);
  M.Psi   = as_matrix(matrices["psi"]);
  M.d     = as_matrix(matrices["thetaDelta"]);
  M.e     = as_matrix(matrices["thetaEpsilon"]);

  M.deriv_dim = 0;
  return M;
}

CompleteData build_complete_data(const Rcpp::List& P,
                                 const Rcpp::List& colidxR,
                                 const arma::uvec& n,
                                 const arma::uvec& d,
                                 const int npatterns) {
  CompleteData data;

  data.V = as_matrix(P["V"]);
  data.npatterns = npatterns;
  data.n = to_eigen(n);
  data.d = to_eigen(d);

  {
    Rcpp::List tgammaR = P["tgamma"];
    data.tgamma.resize(tgammaR.size());
    for (R_xlen_t j = 0; j < tgammaR.size(); ++j) {
      Rcpp::NumericVector vec = tgammaR[j];
      data.tgamma[j] = std::vector<double>(vec.begin(), vec.end());
    }
  }

  {
    Rcpp::List meanR = P["mean"];
    data.mean.resize(meanR.size());
    for (R_xlen_t j = 0; j < meanR.size(); ++j) {
      Rcpp::List inner = meanR[j];
      data.mean[j].reserve(inner.size());
      for (R_xlen_t i = 0; i < inner.size(); ++i) {
        data.mean[j].push_back(as_vector(inner[i]));
      }
    }
  }

  {
    Rcpp::List covR = P["cov"];
    data.cov.resize(covR.size());
    for (R_xlen_t j = 0; j < covR.size(); ++j) {
      Rcpp::List inner = covR[j];
      data.cov[j].reserve(inner.size());
      for (R_xlen_t i = 0; i < inner.size(); ++i) {
        data.cov[j].push_back(as_matrix(inner[i]));
      }
    }
  }

  data.colidx.reserve(colidxR.size());
  for (R_xlen_t i = 0; i < colidxR.size(); ++i) {
    Rcpp::IntegerVector idx = colidxR[i];
    Eigen::VectorXi cols(idx.size());
    for (R_xlen_t j = 0; j < idx.size(); ++j) cols[j] = idx[j];
    data.colidx.push_back(cols);
  }

  return data;
}

ObservedData build_observed_data(const Rcpp::List& dataR,
                                 const Rcpp::List& colidxR,
                                 const Rcpp::List& P,
                                 const arma::uvec& n,
                                 const int npatterns) {
  ObservedData data;
  data.V = as_matrix(P["V"]);
  data.w = as_vector(P["w"]);
  data.n = to_eigen(n);
  data.npatterns = npatterns;

  data.data.reserve(dataR.size());
  for (R_xlen_t i = 0; i < dataR.size(); ++i) {
    data.data.push_back(as_matrix(dataR[i]));
  }

  data.colidx.reserve(colidxR.size());
  for (R_xlen_t i = 0; i < colidxR.size(); ++i) {
    Rcpp::IntegerVector idx = colidxR[i];
    Eigen::VectorXi cols(idx.size());
    for (R_xlen_t j = 0; j < idx.size(); ++j) cols[j] = idx[j];
    data.colidx.push_back(cols);
  }

  return data;
}

template<typename Scalar>
Scalar& model_param(Model<Scalar>& M, std::size_t blk,
                    std::size_t r, std::size_t c) {
  switch (blk) {
    case 0 : return M.lX   (r, c);
    case 1 : return M.lY   (r, c);
    case 2 : return M.tX   (r);
    case 3 : return M.tY   (r);
    case 4 : return M.d    (r, c);
    case 5 : return M.e    (r, c);
    case 6 : return M.A    (r, c);
    case 7 : return M.Psi  (r, c);
    case 8 : return M.a    (r);
    case 9 : return M.beta0(r);
    case 10: return M.Gx   (r, c);
    case 11: return M.Ge   (r, c);
    case 12: return M.Oxx  (r, c);
    case 13: return M.Oex  (r, c);
    default: Rcpp::stop("unknown block id");
  }
}

template<typename Scalar>
Scalar model_param(const Model<Scalar>& M, std::size_t blk,
                   std::size_t r, std::size_t c) {
  return model_param(const_cast<Model<Scalar>&>(M), blk, r, c);
}

template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
get_params_eigen(const Model<Scalar>& M,
                 const arma::uvec& block,
                 const arma::uvec& row,
                 const arma::uvec& col) {
  const std::size_t p = block.n_elem;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> pars(p);
  for (std::size_t k = 0; k < p; ++k) {
    pars(k) = model_param(M, block[k], row[k], col[k]);
  }
  return pars;
}

template<typename Scalar, typename ParamVec>
void set_params_eigen(Model<Scalar>& M,
                      const arma::uvec& block,
                      const arma::uvec& row,
                      const arma::uvec& col,
                      const arma::uvec& symmetric,
                      const ParamVec& vals) {
  const std::size_t p = block.n_elem;
  for (std::size_t k = 0; k < p; ++k) {
    model_param(M, block[k], row[k], col[k]) = vals(k);
    if (symmetric[k] && row[k] != col[k]) {
      model_param(M, block[k], col[k], row[k]) = vals(k);
    }
  }
}

double second_derivative_complete(const Model<double>& base_model,
                                  const CompleteData& data,
                                  const arma::uvec& block,
                                  const arma::uvec& row,
                                  const arma::uvec& col,
                                  const arma::uvec& symmetric,
                                  const Eigen::Matrix<double, Eigen::Dynamic, 1>& theta0,
                                  const Eigen::VectorXd& direction) {
  const std::size_t p = block.n_elem;
  Model<modsem::lms::Dual2> model_dual(base_model);
  set_model_deriv_dim(model_dual, static_cast<Eigen::Index>(p));
  Eigen::Matrix<modsem::lms::Dual2, Eigen::Dynamic, 1> theta(p);
  for (std::size_t i = 0; i < p; ++i)
    theta(i) = modsem::lms::Dual2(theta0(i), direction(i), 0.0);

  set_params_eigen(model_dual, block, row, col, symmetric, theta);

  const modsem::lms::Dual2 ll = modsem::lms::complete_loglik(model_dual, data);
  return ll.d2;
}

double second_derivative_observed(const Model<double>& base_model,
                                  const ObservedData& data,
                                  const arma::uvec& block,
                                  const arma::uvec& row,
                                  const arma::uvec& col,
                                  const arma::uvec& symmetric,
                                  const Eigen::Matrix<double, Eigen::Dynamic, 1>& theta0,
                                  const Eigen::VectorXd& direction) {
  const std::size_t p = block.n_elem;
  Model<modsem::lms::Dual2> model_dual(base_model);
  set_model_deriv_dim(model_dual, static_cast<Eigen::Index>(p));
  Eigen::Matrix<modsem::lms::Dual2, Eigen::Dynamic, 1> theta(p);
  for (std::size_t i = 0; i < p; ++i)
    theta(i) = modsem::lms::Dual2(theta0(i), direction(i), 0.0);

  set_params_eigen(model_dual, block, row, col, symmetric, theta);

  const modsem::lms::Dual2 ll = modsem::lms::observed_loglik(model_dual, data);
  return ll.d2;
}

} // anonymous namespace


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
  (void)eps;
  (void)ncores;

  const Model<double> base_model = build_model_eigen(modelR);
  const CompleteData data = build_complete_data(P, colidxR, n, d, npatterns);

  const std::size_t p = block.n_elem;
  using Der1 = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using AD1  = Eigen::AutoDiffScalar<Der1>;

  Eigen::Matrix<double, Eigen::Dynamic, 1> base =
    get_params_eigen(base_model, block, row, col);

  Eigen::Matrix<AD1, Eigen::Dynamic, 1> theta(p);
  for (std::size_t i = 0; i < p; ++i) {
    theta(i).value() = base(i);
    theta(i).derivatives() = Der1::Zero(p);
    theta(i).derivatives()(i) = 1.0;
  }

  Model<AD1> model_ad(base_model);
  set_model_deriv_dim(model_ad, static_cast<Eigen::Index>(p));
  set_params_eigen(model_ad, block, row, col, symmetric, theta);

  const AD1 ll = modsem::lms::complete_loglik(model_ad, data);

  arma::vec grad(p);
  const Der1& g = ll.derivatives();
  for (std::size_t i = 0; i < p; ++i) grad[i] = g(i);
  return grad;
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
  (void)eps;
  (void)ncores;

  const Model<double> base_model = build_model_eigen(modelR);
  const ObservedData data = build_observed_data(dataR, colidxR, P, n, npatterns);

  const std::size_t p = block.n_elem;
  using Der1 = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using AD1  = Eigen::AutoDiffScalar<Der1>;

  Eigen::Matrix<double, Eigen::Dynamic, 1> base =
    get_params_eigen(base_model, block, row, col);

  Eigen::Matrix<AD1, Eigen::Dynamic, 1> theta(p);
  for (std::size_t i = 0; i < p; ++i) {
    theta(i).value() = base(i);
    theta(i).derivatives() = Der1::Zero(p);
    theta(i).derivatives()(i) = 1.0;
  }

  Model<AD1> model_ad(base_model);
  set_model_deriv_dim(model_ad, static_cast<Eigen::Index>(p));
  set_params_eigen(model_ad, block, row, col, symmetric, theta);

  const AD1 ll = modsem::lms::observed_loglik(model_ad, data);

  arma::vec grad(p);
  const Der1& g = ll.derivatives();
  for (std::size_t i = 0; i < p; ++i) grad[i] = g(i);
  return grad;
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
    set_params(Mc, block, row, col, symmetric, base + disp[k] % incr);
    y[k] = fun(Mc);
  }

  // Restore baseline
  set_params(M, block, row, col, symmetric, base);

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
  auto pair_index = [p](std::size_t i, std::size_t j) -> std::size_t {
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
        std::size_t k = pair_index(i,j);
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
    set_params(Mc, block, row, col, symmetric, base + disp[k] % incr);
    y[k] = fun(Mc);
  }
  set_params(M, block, row, col, symmetric, base);

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
        std::size_t k = pair_index(i,j);
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
  const arma::vec base = get_params(M, block, row, col);
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
  (void)relStep;
  (void)minAbs;
  (void)ncores;

  const Model<double> base_model = build_model_eigen(modelR);
  const ObservedData data = build_observed_data(dataR, colidxR, P, n, npatterns);

  const std::size_t p = block.n_elem;
  Eigen::Matrix<double, Eigen::Dynamic, 1> base =
    get_params_eigen(base_model, block, row, col);

  // Compute gradient via first-order AD
  arma::vec grad(p);
  {
    using Der1 = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using AD1  = Eigen::AutoDiffScalar<Der1>;
    Eigen::Matrix<AD1, Eigen::Dynamic, 1> theta(p);
    for (std::size_t i = 0; i < p; ++i) {
      theta(i).value() = base(i);
      theta(i).derivatives() = Der1::Zero(p);
      theta(i).derivatives()(i) = 1.0;
    }
    Model<AD1> model_ad(base_model);
    set_model_deriv_dim(model_ad, static_cast<Eigen::Index>(p));
    set_params_eigen(model_ad, block, row, col, symmetric, theta);
    const AD1 ll = modsem::lms::observed_loglik(model_ad, data);
    const Der1& g = ll.derivatives();
    for (std::size_t i = 0; i < p; ++i) grad[i] = g(i);
  }

  arma::mat Hess(p, p, arma::fill::zeros);
  Eigen::VectorXd direction = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd diag(p);

  for (std::size_t i = 0; i < p; ++i) {
    direction.setZero();
    direction(i) = 1.0;
    const double d2 = second_derivative_observed(base_model, data,
      block, row, col, symmetric, base, direction);
    Hess(i, i) = d2;
    diag(i) = d2;
  }

  for (std::size_t i = 0; i < p; ++i) {
    for (std::size_t j = i + 1; j < p; ++j) {
      direction.setZero();
      direction(i) = 1.0;
      direction(j) = 1.0;
      const double d2 = second_derivative_observed(base_model, data,
        block, row, col, symmetric, base, direction);
      const double off = 0.5 * (d2 - diag(i) - diag(j));
      Hess(i, j) = off;
      Hess(j, i) = off;
    }
  }

  const double value = modsem::lms::observed_loglik(base_model, data);

  return Rcpp::List::create(
    Rcpp::Named("mean")     = value,
    Rcpp::Named("gradient") = grad,
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
  (void)relStep;
  (void)minAbs;
  (void)ncores;

  const Model<double> base_model = build_model_eigen(modelR);
  const CompleteData data = build_complete_data(P, colidxR, n, d, npatterns);

  const std::size_t p = block.n_elem;
  Eigen::Matrix<double, Eigen::Dynamic, 1> base =
    get_params_eigen(base_model, block, row, col);

  arma::vec grad(p);
  {
    using Der1 = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using AD1  = Eigen::AutoDiffScalar<Der1>;
    Eigen::Matrix<AD1, Eigen::Dynamic, 1> theta(p);
    for (std::size_t i = 0; i < p; ++i) {
      theta(i).value() = base(i);
      theta(i).derivatives() = Der1::Zero(p);
      theta(i).derivatives()(i) = 1.0;
    }
    Model<AD1> model_ad(base_model);
    set_model_deriv_dim(model_ad, static_cast<Eigen::Index>(p));
    set_params_eigen(model_ad, block, row, col, symmetric, theta);
    const AD1 ll = modsem::lms::complete_loglik(model_ad, data);
    const Der1& g = ll.derivatives();
    for (std::size_t i = 0; i < p; ++i) grad[i] = g(i);
  }

  arma::mat Hess(p, p, arma::fill::zeros);
  Eigen::VectorXd direction = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd diag(p);

  for (std::size_t i = 0; i < p; ++i) {
    direction.setZero();
    direction(i) = 1.0;
    const double d2 = second_derivative_complete(base_model, data,
      block, row, col, symmetric, base, direction);
    Hess(i, i) = d2;
    diag(i) = d2;
  }

  for (std::size_t i = 0; i < p; ++i) {
    for (std::size_t j = i + 1; j < p; ++j) {
      direction.setZero();
      direction(i) = 1.0;
      direction(j) = 1.0;
      const double d2 = second_derivative_complete(base_model, data,
        block, row, col, symmetric, base, direction);
      const double off = 0.5 * (d2 - diag(i) - diag(j));
      Hess(i, j) = off;
      Hess(j, i) = off;
    }
  }

  const double value = modsem::lms::complete_loglik(base_model, data);

  return Rcpp::List::create(
    Rcpp::Named("mean")     = value,
    Rcpp::Named("gradient") = grad,
    Rcpp::Named("Hessian")  = Hess
  );
}
