#include <RcppArmadillo.h>
#include "lms.h"
#include "mvnorm.h"
#include <float.h>

// [[Rcpp::depends(RcppArmadillo)]]


// Deprecated will remove soon...
// [[Rcpp::export]]
arma::vec muLmsCpp(Rcpp::List model, arma::vec z) {
  Rcpp::List matrices = model["matrices"];
  Rcpp::List info = model["info"];
  Rcpp::List quad = model["quad"];
  int numXis = Rcpp::as<int>(info["numXis"]);
  int k = Rcpp::as<int>(quad["k"]);
  arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  arma::mat Ie = Rcpp::as<arma::mat>(matrices["Ieta"]);
  arma::mat lY = Rcpp::as<arma::mat>(matrices["lambdaY"]);
  arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  arma::mat tY = Rcpp::as<arma::mat>(matrices["tauY"]);
  arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else zVec = arma::zeros<arma::vec>(numXis);
  arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);
  }
  arma::vec muX = tX + lX * (beta0 + A * zVec);
  arma::vec muY = tY + 
    lY * (Binv * (a + 
          Gx * (beta0 + A * zVec) + 
          kronZ.t() * Oxx * (beta0 + A * zVec)));
  return arma::join_cols(muX, muY);
}


// Deprecated will remove soon...
// [[Rcpp::export]]
arma::mat sigmaLmsCpp(Rcpp::List model, arma::vec z) {
  Rcpp::List matrices = model["matrices"];
  Rcpp::List info = model["info"];
  Rcpp::List quad = model["quad"];
  int numXis = Rcpp::as<int>(info["numXis"]);
  int k = Rcpp::as<int>(quad["k"]);
  arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  arma::mat Ie = Rcpp::as<arma::mat>(matrices["Ieta"]);
  arma::mat lY = Rcpp::as<arma::mat>(matrices["lambdaY"]);
  arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  arma::mat tY = Rcpp::as<arma::mat>(matrices["tauY"]);
  arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);
  arma::mat Psi = Rcpp::as<arma::mat>(matrices["psi"]); 
  arma::mat d = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  arma::mat e = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else zVec = arma::zeros<arma::vec>(numXis);
  arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);
  }

  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));
  arma::mat Sxx = lX * A * Oi * A.t() * lX.t() + d;
  arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
  arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
  arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() + 
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
  else return arma::zeros<arma::vec>(numXis);
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
    arma::vec zVec = make_zvec(k, numXis, z);
    arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

    arma::mat Binv;
    if (Ie.n_cols == 1) Binv = arma::mat(Ie);
    else Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    arma::vec muX = tX + lX * (beta0 + A * zVec);
    arma::vec muY = tY + 
      lY * (Binv * (a + 
            Gx * (beta0 + A * zVec) + 
            kronZ.t() * Oxx * (beta0 + A * zVec)));
    return arma::join_cols(muX, muY);
  }


  arma::mat Sigma(const arma::vec& z) const {
    arma::vec zVec = make_zvec(k, numXis, z);

    arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

    arma::mat Binv;
    if (Ie.n_cols == 1) Binv = arma::mat(Ie);
    else Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    arma::mat Oi = make_Oi(k, numXis);
    arma::mat Sxx = lX * A * Oi * A.t() * lX.t() + d;
    arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
    arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
    arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() + 
      lY * (Binv * Psi * Binv.t()) * lY.t() + e;

    return arma::join_cols(
        arma::join_rows(Sxx, Sxy), 
        arma::join_rows(Sxy.t(), Syy)
        );
  }
};



inline std::vector<arma::vec> as_vec_of_vec(const Rcpp::List& L) {
    const std::size_t J = L.size();
    std::vector<arma::vec> out;
    out.reserve(J);
    for (std::size_t j = 0; j < J; ++j)
        out.emplace_back( Rcpp::as<arma::vec>(L[j]) );
    return out;
}


inline std::vector<arma::mat> as_vec_of_mat(const Rcpp::List& L) {
    const std::size_t J = L.size();
    std::vector<arma::mat> out;
    out.reserve(J);
    for (std::size_t j = 0; j < J; ++j)
        out.emplace_back( Rcpp::as<arma::mat>(L[j]) );
    return out;
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
    default: Rcpp::stop("unknown block id");
  }
}


template< class F >
arma::vec gradientFD(LMSModel&         M,       
                     F&&               logLik,  
                     const arma::uvec& block,
                     const arma::uvec& row,
                     const arma::uvec& col,
                     const double      eps = 1e-6) {
  const std::size_t p = block.n_elem;
  arma::vec grad(p);

  const double f0 = logLik(M);

  for (std::size_t k = 0; k < p; ++k) {

    double& ti   = lms_param(M, block[k], row[k], col[k]);
    const double old = ti;

    ti += eps;                  
    const double f1 = logLik(M);

    grad[k] = (f1 - f0) / eps;

    ti = old;
  }

  return grad;
}


inline double
completeLogLikFromModel(const LMSModel&                 M,
                        const arma::mat&                V,          // J × k
                        const arma::vec&                TGamma,     // length J
                        const std::vector<arma::vec>&   Mean,       // length J
                        const std::vector<arma::mat>&   Cov,        // length J
                        const int                       n,
                        const int                       d) {
    const std::size_t J = V.n_rows;
    double ll = 0.0;

    for (std::size_t j = 0; j < J; ++j) {

        const double tg = TGamma[j];
        if (tg <= DBL_MIN) continue;

        const arma::vec& z  = V.row(j).t();   // view – no copy
        const arma::vec& nu = Mean[j];
        const arma::mat& S  = Cov [j];

        const arma::vec mu  = M.mu   (z);
        const arma::mat Sig = M.Sigma(z);

        ll += totalDmvnWeightedCpp(mu, Sig, nu, S, tg, n, d);
    }
    return ll;
}


// [[Rcpp::export]]
arma::vec gradLogLikLmsCpp(const Rcpp::List& modelR,
                           const Rcpp::List& P,
                           const arma::uvec& block,
                           const arma::uvec& row,
                           const arma::uvec& col,
                           double            eps = 1e-6) {
    LMSModel M(modelR);

    const arma::mat  V       = Rcpp::as<arma::mat>(P["V"]);
    const arma::vec  tgamma  = Rcpp::as<arma::vec>(P["tgamma"]);
    const auto       Mean    = as_vec_of_vec(P["mean"]);
    const auto       Cov     = as_vec_of_mat(P["cov"]);

    const Rcpp::List info   = modelR["info"];
    const int n             = Rcpp::as<int>(info["N"]);
    const int d             = Rcpp::as<int>(info["ncol"]);

    auto comp_ll = [&](LMSModel& mod) -> double {
        return completeLogLikFromModel(mod, V, tgamma, Mean, Cov, n, d);
    };

    return gradientFD(M, comp_ll, block, row, col, eps);
}


// [[Rcpp::export]]
double completeLogLikLmsCpp(Rcpp::List modelR, Rcpp::List P, Rcpp::List quad) {
  const LMSModel model = LMSModel(modelR);
  const arma::mat Z = Rcpp::as<arma::mat>(quad["n"]).t(); // transpose so we can use column order vectors

  const arma::mat  V       = Rcpp::as<arma::mat>(P["V"]);
  const arma::vec  tgamma  = Rcpp::as<arma::vec>(P["tgamma"]);
  const auto       Mean    = as_vec_of_vec(P["mean"]);
  const auto       Cov     = as_vec_of_mat(P["cov"]);

  const Rcpp::List info   = modelR["info"];
  const int n             = Rcpp::as<int>(info["N"]);
  const int d             = Rcpp::as<int>(info["ncol"]);

  return completeLogLikFromModel(model, V, tgamma, Mean, Cov, n, d);
}


inline double observedLogLikFromModel(const LMSModel&            M,
                                      const arma::mat&           V,
                                      const arma::vec&           w, 
                                      const arma::mat&        data,
                                      const int ncores = 1) {
    const std::size_t n = V.n_rows;

    arma::vec density = arma::zeros<arma::vec>(data.n_rows);
    for (std::size_t i = 0; i < n; ++i) {
      if (w[i] <= DBL_MIN) continue;

      const arma::vec z   = V.row(i).t();
      const arma::vec mu  = M.mu   (z);
      const arma::mat Sig = M.Sigma(z);

      density += dmvnfast(data, mu, Sig, false, ncores, false) * w[i];
    }

    return arma::sum(arma::log(density));
}


// [[Rcpp::export]]
arma::vec gradObsLogLikLmsCpp(const Rcpp::List& modelR,
                              const arma::mat&  data,
                              const Rcpp::List& P,
                              const arma::uvec& block,
                              const arma::uvec& row,
                              const arma::uvec& col,
                              double            eps   = 1e-6,
                              int               ncores= 1) {
  LMSModel M(modelR);

  const arma::mat V = Rcpp::as<arma::mat>(P["V"]);
  const arma::vec w = Rcpp::as<arma::vec>(P["w"]);

  auto obs_ll = [&](LMSModel& mod) -> double {
    return observedLogLikFromModel(mod, V, w, data, ncores);
  };

  return gradientFD(M, obs_ll, block, row, col, eps);
}


// [[Rcpp::export]]
double observedLogLikLmsCpp(Rcpp::List modelR, arma::mat data, Rcpp::List P, const int ncores = 1) {
  const LMSModel M = LMSModel(modelR);

  const arma::mat V       = Rcpp::as<arma::mat>(P["V"]);
  const arma::vec w       = Rcpp::as<arma::vec>(P["w"]);

  return observedLogLikFromModel(M, V, w, data, ncores);
}
