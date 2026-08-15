#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

#include "lms_model.h"
// [[Rcpp::depends(RcppArmadillo)]]

// The PML kernel: implied moments, bivariate normal CDF, rectangles, and the
// composite log-likelihood.
//
// Conditional on the nonlinear latents the LMS model is linear-Gaussian, so
// `LMSModel::muSigma(z)` gives the mean and covariance of the underlying
// variables at a node and a pair probability is a rectangle under the
// corresponding bivariate marginal. Going through LMSModel rather than
// reimplementing those moments means PML inherits every fix the LMS kernel
// gets, composites included -- and the R version this replaced HAD already
// drifted: it conditioned on xi values while being fed standard normal
// innovations, ignored omegaEtaXi, and only ever worked for a single eta.
//
// Everything here is called from R/est_pml.R once per objective evaluation.
// The data enters only as contingency tables, so N is not in the cost.


// ---------------------------------------------------------------------------
// Bivariate normal CDF, after Genz (2004) -- the same algorithm mvtnorm uses.
//
// PML lives or dies on this: a fit evaluates it millions of times, so calling
// out to R was never viable. `bvnUpper` returns P(X > h, Y > k) and switches
// between three Gauss-Legendre rules by |rho|, with a separate near-singular
// expansion above 0.925 where the plain quadrature loses accuracy.

static const double kTwoPi = 6.283185307179586;

static const double kW6[3]  = {0.1713244923791705, 0.3607615730481384,
                               0.4679139345726904};
static const double kX6[3]  = {0.9324695142031522, 0.6612093864662647,
                               0.2386191860831970};
static const double kW12[6] = {0.04717533638651177, 0.1069393259953183,
                               0.1600783285433464, 0.2031674267230659,
                               0.2334925365383547, 0.2491470458134029};
static const double kX12[6] = {0.9815606342467191, 0.9041172563704750,
                               0.7699026741943050, 0.5873179542866171,
                               0.3678314989981802, 0.1252334085114692};
static const double kW20[10] = {0.01761400713915212, 0.04060142980038694,
                                0.06267204833410906, 0.08327674157670475,
                                0.1019301198172404,  0.1181945319615184,
                                0.1316886384491766,  0.1420961093183821,
                                0.1491729864726037,  0.1527533871307259};
static const double kX20[10] = {0.9931285991850949, 0.9639719272779138,
                                0.9122344282513259, 0.8391169718222188,
                                0.7463319064601508, 0.6360536807265150,
                                0.5108670019508271, 0.3737060887154196,
                                0.2277858511416451, 0.07652652113349733};

static inline double normalCdf(double z) {
  return R::pnorm(z, 0.0, 1.0, 1, 0);
}


static double bvnUpper(double dh, double dk, double r) {
  double h = dh, k = dk, hk = h * k, bvn = 0.0;

  int lg;
  const double *w, *x;
  if (std::fabs(r) < 0.3)       { lg = 3;  w = kW6;  x = kX6;  }
  else if (std::fabs(r) < 0.75) { lg = 6;  w = kW12; x = kX12; }
  else                          { lg = 10; w = kW20; x = kX20; }

  if (std::fabs(r) < 0.925) {
    const double hs  = (h * h + k * k) / 2.0;
    const double asr = std::asin(r);
    for (int i = 0; i < lg; ++i) {
      for (int is = -1; is <= 1; is += 2) {
        const double sn = std::sin(asr * (is * x[i] + 1.0) / 2.0);
        bvn += w[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
      }
    }
    bvn = bvn * asr / (2.0 * kTwoPi) + normalCdf(-h) * normalCdf(-k);

  } else {
    if (r < 0) { k = -k; hk = -hk; }

    if (std::fabs(r) < 1) {
      const double as = (1.0 - r) * (1.0 + r);
      double a        = std::sqrt(as);
      const double bs = (h - k) * (h - k);
      const double c  = (4.0 - hk) / 8.0;
      const double d  = (12.0 - hk) / 16.0;

      double asr = -(bs / as + hk) / 2.0;
      if (asr > -100.0)
        bvn = a * std::exp(asr) *
              (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 +
               c * d * as * as / 5.0);

      if (hk > -100.0) {
        const double b = std::sqrt(bs);
        bvn -= std::exp(-hk / 2.0) * std::sqrt(kTwoPi) * normalCdf(-b / a) * b *
               (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
      }

      a /= 2.0;
      for (int i = 0; i < lg; ++i) {
        for (int is = -1; is <= 1; is += 2) {
          double xs = a * (is * x[i] + 1.0);
          xs = xs * xs;
          const double rs = std::sqrt(1.0 - xs);
          const double e  = -(bs / xs + hk) / 2.0;
          if (e > -100.0)
            bvn += a * w[i] * std::exp(e) *
                   (std::exp(-hk * xs / (2.0 * (1.0 + rs) * (1.0 + rs))) / rs -
                    (1.0 + c * xs * (1.0 + d * xs)));
        }
      }
      bvn = -bvn / kTwoPi;
    }

    if (r > 0) {
      bvn += normalCdf(-std::max(h, k));
    } else {
      bvn = -bvn;
      if (k > h) bvn += normalCdf(k) - normalCdf(h);
    }
  }

  if (bvn < 0.0) bvn = 0.0;
  if (bvn > 1.0) bvn = 1.0;
  return bvn;
}


// P(X <= a, Y <= b), infinite arguments included -- every rectangle's outer row
// and column needs them.
static inline double bvnLower(double a, double b, double r) {
  if (a == R_NegInf || b == R_NegInf) return 0.0;
  if (a == R_PosInf) return b == R_PosInf ? 1.0 : normalCdf(b);
  if (b == R_PosInf) return normalCdf(a);
  return bvnUpper(-a, -b, r);
}


// Exposed for testing against pbivnorm; not used by the objective.
// [[Rcpp::export]]
arma::vec pmlBivariateCpp(const arma::vec& a, const arma::vec& b, double rho) {
  arma::vec out(a.n_elem);
  for (arma::uword i = 0; i < a.n_elem; ++i) out[i] = bvnLower(a[i], b[i], rho);
  return out;
}


// ---------------------------------------------------------------------------
// Cell probabilities for one pair.
//
// Every cell of the table reads from the same grid of corners, so the corners
// are evaluated once -- (Kj+1)(Kk+1) CDF calls instead of 4*Kj*Kk -- and the
// cells come out as second differences. rho is fixed across the whole grid,
// which is what makes bvnUpper's rule selection amortise.

static void pairProbabilities(const arma::vec& mu, const arma::mat& sigma,
                              const std::vector<arma::vec>& thresholds,
                              const arma::uvec& rows,
                              arma::uword j, arma::uword k,
                              double weight, bool accumulate,
                              arma::mat& out) {
  const arma::uword rj = rows[j], rk = rows[k];
  const double sj  = std::sqrt(sigma(rj, rj));
  const double sk  = std::sqrt(sigma(rk, rk));
  const double rho = sigma(rj, rk) / (sj * sk);

  const arma::vec& tj = thresholds[j];
  const arma::vec& tk = thresholds[k];
  const arma::uword nj = tj.n_elem, nk = tk.n_elem;

  std::vector<double> a(nj + 2), b(nk + 2);
  a[0] = R_NegInf; a[nj + 1] = R_PosInf;
  b[0] = R_NegInf; b[nk + 1] = R_PosInf;
  for (arma::uword c = 0; c < nj; ++c) a[c + 1] = (tj[c] - mu[rj]) / sj;
  for (arma::uword d = 0; d < nk; ++d) b[d + 1] = (tk[d] - mu[rk]) / sk;

  arma::mat corner(nj + 2, nk + 2);
  for (arma::uword c = 0; c < nj + 2; ++c)
    for (arma::uword d = 0; d < nk + 2; ++d)
      corner(c, d) = bvnLower(a[c], b[d], rho);

  for (arma::uword c = 0; c < nj + 1; ++c) {
    for (arma::uword d = 0; d < nk + 1; ++d) {
      const double cell = corner(c + 1, d + 1) - corner(c, d + 1) -
                          corner(c + 1, d)     + corner(c, d);
      if (accumulate) out(c, d) += weight * cell;
      else            out(c, d)  = weight * cell;
    }
  }
}


// Cell probabilities for every pair.
//
// `hoisted` pairs are evaluated once from the unconditional moments and
// `integrated` pairs inside the node loop; R decides which is which, from the
// model's structure, before the fit starts. All indices are 0-based, and
// `rows` maps the ordered indicators onto rows of mu/Sigma.
static std::vector<arma::mat> pmlProbabilities(
    const Rcpp::List& modFilled, const arma::mat& nodes,
    const arma::vec& weights, const std::vector<arma::vec>& thresholds,
    const arma::uvec& rows, const arma::umat& pairs,
    const arma::uvec& hoisted, const arma::uvec& integrated) {
  LMSModel model(modFilled);

  std::vector<arma::mat> probability(pairs.n_rows);
  for (arma::uword p = 0; p < pairs.n_rows; ++p)
    probability[p] = arma::zeros<arma::mat>(thresholds[pairs(p, 0)].n_elem + 1,
                                            thresholds[pairs(p, 1)].n_elem + 1);

  if (hoisted.n_elem) {
    LMSModel marginal = model;
    marginal.k = 0;
    marginal.updateCache();
    const std::pair<arma::vec, arma::mat> ms = marginal.muSigma(arma::vec());
    for (arma::uword i = 0; i < hoisted.n_elem; ++i) {
      const arma::uword p = hoisted[i];
      pairProbabilities(ms.first, ms.second, thresholds, rows,
                        pairs(p, 0), pairs(p, 1), 1.0, false, probability[p]);
    }
  }

  for (arma::uword q = 0; q < nodes.n_rows && integrated.n_elem; ++q) {
    const std::pair<arma::vec, arma::mat> ms = model.muSigma(nodes.row(q).t());
    for (arma::uword i = 0; i < integrated.n_elem; ++i) {
      const arma::uword p = integrated[i];
      pairProbabilities(ms.first, ms.second, thresholds, rows,
                        pairs(p, 0), pairs(p, 1), weights[q], true,
                        probability[p]);
    }
  }

  return probability;
}


static std::vector<arma::vec> asThresholds(const Rcpp::List& thresholdList) {
  std::vector<arma::vec> out(thresholdList.size());
  for (int i = 0; i < thresholdList.size(); ++i)
    out[i] = Rcpp::as<arma::vec>(thresholdList[i]);
  return out;
}


// The pairwise composite log-likelihood.
// [[Rcpp::export]]
double pmlObjectiveCpp(const Rcpp::List& modFilled,
                       const arma::mat& nodes,
                       const arma::vec& weights,
                       const Rcpp::List& thresholdList,
                       const arma::uvec& rows,
                       const arma::umat& pairs,
                       const Rcpp::List& countList,
                       const arma::uvec& hoisted,
                       const arma::uvec& integrated) {
  const std::vector<arma::mat> probability = pmlProbabilities(
    modFilled, nodes, weights, asThresholds(thresholdList), rows, pairs,
    hoisted, integrated);

  double total = 0.0;
  for (arma::uword p = 0; p < pairs.n_rows; ++p) {
    const arma::mat n = Rcpp::as<arma::mat>(countList[p]);
    const arma::mat& probs = probability[p];
    for (arma::uword c = 0; c < n.n_rows; ++c)
      for (arma::uword d = 0; d < n.n_cols; ++d)
        if (n(c, d) > 0)
          total += n(c, d) * std::log(std::max(probs(c, d), 1e-12));
  }

  return total;
}


// The same probabilities, handed back for inspection. This is what the tests
// check against Monte Carlo and against lavaan, so they exercise the kernel the
// fit actually runs rather than a second implementation of it.
// [[Rcpp::export]]
Rcpp::List pmlProbabilitiesCpp(const Rcpp::List& modFilled,
                               const arma::mat& nodes,
                               const arma::vec& weights,
                               const Rcpp::List& thresholdList,
                               const arma::uvec& rows,
                               const arma::umat& pairs,
                               const arma::uvec& hoisted,
                               const arma::uvec& integrated) {
  const std::vector<arma::mat> probability = pmlProbabilities(
    modFilled, nodes, weights, asThresholds(thresholdList), rows, pairs,
    hoisted, integrated);
  return Rcpp::wrap(probability);
}
