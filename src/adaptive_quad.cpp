// -----------------------------------------------------------------------------
//  lms_aghq.cpp – Rcpp/Armadillo backend for adaptive‑quadrature LMS
//  All Armadillo types fully qualified (no namespace aliases).
// -----------------------------------------------------------------------------
#include <RcppArmadillo.h>
#include <numeric>
#include <cfloat>
#include <limits>        // for quiet_NaN()
// [[Rcpp::depends(RcppArmadillo)]]

static const double log2pi = std::log(2.0 * M_PI);
static const double kNaN   = std::numeric_limits<double>::quiet_NaN();

inline double dmvnorm_fast(const arma::vec& x,
                           const arma::vec& mu,
                           const arma::mat& S) {

    // 1. Try a *fast* symmetric-positive-definite test.
    //    If it fails, short-circuit to NaN.  No exception is thrown.
    arma::mat L;
    if (!arma::chol(L, S, "lower"))                 // false ⇒ S not SPD
        return kNaN;

    // 2. Everything is fine – proceed with the Cholesky-based evaluation.
    const arma::uword d = x.n_elem;

    const double log_det = 2.0 * arma::sum(log(L.diag()));  // log|S|

    arma::vec diff = x - mu;

    // Solve L * y = diff  →  y = L⁻¹ diff   (no explicit inverse)
    arma::vec y = arma::solve(arma::trimatl(L), diff);

    // Quadratic form  (x−μ)ᵀ S⁻¹ (x−μ)  =  yᵀ y
    const double quad = arma::dot(y, y);

    // 3. Catch any NaNs produced by sub-operations (paranoia guard)
    if (!std::isfinite(log_det) || !std::isfinite(quad))
        return kNaN;

    return std::exp(-0.5 * (d * log2pi + log_det + quad));
}


class LmsModel {
public:
  arma::mat A, omegaXiXi, omegaEtaXi, Ieta,
            lambdaY, lambdaX, tauY, tauX,
            gammaXi, gammaEta, alpha, psi,
            thetaDelta, thetaEpsilon;
  arma::uword k, numXis;

  explicit LmsModel(const Rcpp::List& mod) {
    Rcpp::List M = mod["matrices"];
#define GRAB(name) name = Rcpp::as<arma::mat>(M[#name]);
    GRAB(A) GRAB(omegaXiXi) GRAB(omegaEtaXi) GRAB(Ieta)
    GRAB(lambdaY) GRAB(lambdaX) GRAB(tauY) GRAB(tauX)
    GRAB(gammaXi) GRAB(gammaEta) GRAB(alpha) GRAB(psi)
    GRAB(thetaDelta) GRAB(thetaEpsilon)
#undef GRAB
    Rcpp::List q = mod["quad"], info = mod["info"];
    k      = q.containsElementNamed("k")      ? static_cast<arma::uword>(Rcpp::as<int>(q["k"]))       : 1;
    numXis = info.containsElementNamed("numXis") ? static_cast<arma::uword>(Rcpp::as<int>(info["numXis"])) : k;
  }

  arma::vec mu(const arma::vec& z1) const {
    arma::vec z(numXis, arma::fill::zeros); z(arma::span(0, k-1)) = z1;
    arma::mat kronZ = arma::kron(Ieta, A * z);
    arma::mat Binv  = (Ieta.n_cols == 1) ? Ieta
      : arma::inv(Ieta - gammaEta - kronZ.t() * omegaEtaXi);

    arma::vec muX = tauX + lambdaX * A * z;
    arma::vec muY = tauY + lambdaY * (Binv * (alpha + gammaXi * A * z + kronZ.t() * omegaXiXi * A * z));
    return arma::join_vert(muX, muY);
  }

  arma::mat sigma(const arma::vec& z1) const {
    arma::vec z(numXis, arma::fill::zeros); z(arma::span(0, k-1)) = z1;
    arma::mat kronZ = arma::kron(Ieta, A * z);
    arma::mat Binv  = (Ieta.n_cols == 1) ? Ieta
      : arma::inv(Ieta - gammaEta - kronZ.t() * omegaEtaXi);

    arma::mat OI = arma::eye(numXis, numXis);
    OI(arma::span(0, k-1), arma::span(0, k-1)).zeros();

    arma::mat Sxx = lambdaX * A * OI * A.t() * lambdaX.t() + thetaDelta;
    arma::mat G   = Binv * (gammaXi * A + kronZ.t() * omegaXiXi * A);
    arma::mat Sxy = lambdaX * A * OI * G.t() * lambdaY.t();
    arma::mat Syy = lambdaY * (G * OI * G.t() + Binv * psi * Binv.t()) * lambdaY.t() + thetaEpsilon;

    arma::mat out(Sxx.n_rows + Syy.n_rows, Sxx.n_rows + Syy.n_rows, arma::fill::zeros);
    out(arma::span(0, Sxx.n_rows-1), arma::span(0, Sxx.n_cols-1)) = Sxx;
    out(arma::span(0, Sxy.n_rows-1), arma::span(Sxx.n_cols, out.n_cols-1)) = Sxy;
    out(arma::span(Sxx.n_rows, out.n_rows-1), arma::span(0, Sxy.n_rows-1)) = Sxy.t();
    out(arma::span(Sxx.n_rows, out.n_rows-1), arma::span(Sxx.n_cols, out.n_cols-1)) = Syy;
    return out;
  }
};


struct QuadCache {
  std::vector<arma::mat> V; 
  std::vector<arma::vec> W;
  std::vector<arma::vec> P; 

  explicit QuadCache(std::size_t N): V(N), W(N), P(N) {}
};


using quadPtr = Rcpp::XPtr<QuadCache>;


// [[Rcpp::export]]
SEXP QuadGrid2XPtr(const Rcpp::List &quadR) {
  const Rcpp::List P = quadR["P"];
  const Rcpp::List V = quadR["V"];
  const Rcpp::List W = quadR["w"];
  
  const int N = P.size();
  QuadCache *cache = new QuadCache(N);

  for (int i = 0; i < N; i++) {
    cache->V[i] = Rcpp::as<arma::mat>(V[i]);
    cache->W[i] = Rcpp::as<arma::vec>(W[i]);
    cache->P[i] = Rcpp::as<arma::vec>(P[i]);
  }

  return Rcpp::XPtr<QuadCache>(cache, true); // let Rcpp garbage collect the memory
}


// [[Rcpp::export]]
double completeLogLikLmsCpp(const Rcpp::List& modelR, const arma::mat& data, const SEXP xpCache) {
  LmsModel mod(modelR);

  const int N = data.n_rows;
  Rcpp::XPtr<QuadCache> cache(xpCache); 

  double ll = 0;

  for (int i = 0; i < N; i++) {
    const arma::vec y_i = data.row(i).t();
    const arma::vec p_i = cache->P[i]; // posterior probs for case i (length m_i)
    const arma::mat z_i = cache->V[i]; // adaptive nodes (m_i × k)
    const int m_i = z_i.n_rows;

    for (int j = 0; j < m_i; j++) {
      const double p_ij = p_i(j);

      if (p_ij <= DBL_MIN) continue;

      const arma::rowvec z_ij = z_i.row(j);
      const arma::vec mu_ij = mod.mu(z_ij);
      const arma::mat sigma_ij = mod.sigma(z_ij); 
 
      const double ll_ij = std::log(dmvnorm_fast(y_i, mu_ij, sigma_ij));

      if (std::isnan(ll_ij)) return kNaN;

      ll += p_ij * ll_ij;
    }
  }

  return ll;
}


// [[Rcpp::export]]
double obsLogLikLmsCpp(const Rcpp::List& modelR, const arma::mat& data,
    const Rcpp::List& quadR) {
  LmsModel mod(modelR);

  const int N = data.n_rows;

  Rcpp::List P = quadR["P"];
  Rcpp::List V = quadR["V"];
  Rcpp::List W = quadR["w"];

  double ll = 0;

  for (int i = 0; i < N; i++) {
    const arma::vec y_i = data.row(i).t();
    const arma::vec w_i = Rcpp::as<arma::vec>(W[i]);           // posterior probs for case i (length m_i)
    const arma::mat z_i = Rcpp::as<arma::mat>(V[i]);             // adaptive nodes (m_i × k)
    const int m_i = z_i.n_rows;

    arma::vec contrib(m_i);

    for (int j = 0; j < m_i; j++) {
      const double w_ij = w_i(j);
      const arma::rowvec z_ij = z_i.row(j);
      const arma::vec mu_ij = mod.mu(z_ij);
      const arma::mat sigma_ij = mod.sigma(z_ij); 

      contrib(j) = w_ij * dmvnorm_fast(y_i, mu_ij, sigma_ij);
    }

    ll += std::log(arma::sum(contrib));
  }

  return ll;
}
