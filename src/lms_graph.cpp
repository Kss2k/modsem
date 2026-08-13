#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct GraphTerm {
  int first;
  int second;
  double coefficient;
};

// One factor of a product term: a latent variable and how many times it
// appears, so `X:X` is a single factor with exponent 2.
struct ProductFactor {
  int variable;
  int exponent;
};

// A product is described once, globally, and referenced by index from every
// equation that uses it, so it is evaluated at most once per node.
struct ProductPlan {
  std::vector<ProductFactor> factors;
  int lastEta; // -1 if the product contains no endogenous variable
};

struct GraphProductTerm {
  int product;
  double coefficient;
};

struct GraphEquationPlan {
  std::vector<GraphTerm> xiLinear;
  std::vector<GraphTerm> etaLinear;
  std::vector<GraphProductTerm> products;
};

struct GraphPlan {
  std::vector<GraphEquationPlan> equations;
  std::vector<ProductPlan> products;
};

struct LoadingTerm {
  int indicator;
  int latent;
  double coefficient;
};

std::vector<LoadingTerm> buildLoadingPlan(const arma::mat& lambda) {
  std::vector<LoadingTerm> plan;
  plan.reserve(arma::accu(lambda != 0.0));
  for (arma::uword indicator = 0; indicator < lambda.n_rows; ++indicator)
    for (arma::uword latent = 0; latent < lambda.n_cols; ++latent)
      if (lambda(indicator, latent) != 0.0)
        plan.push_back({static_cast<int>(indicator), static_cast<int>(latent),
                        lambda(indicator, latent)});
  return plan;
}

arma::mat measurementMeans(const arma::mat& states,
                           const arma::vec& intercepts,
                           const std::vector<LoadingTerm>& plan) {
  arma::mat means(states.n_rows, intercepts.n_elem, arma::fill::zeros);
  means.each_row() += intercepts.t();
  for (const LoadingTerm& term : plan)
    means.col(term.indicator) += term.coefficient * states.col(term.latent);
  return means;
}

// `productDesign` is P by (numXis + numEtas) integer exponents; `omega` is
// numEtas by P coefficients. Together they replace the stacked Kronecker omega
// blocks, which could only express two-way products and forbade eta-by-eta
// terms. See R/lms_graph_products.R for how they are built.
GraphPlan buildGraphPlan(const arma::mat& gammaXi,
                         const arma::mat& gammaEta,
                         const arma::mat& productDesign,
                         const arma::mat& omega,
                         const int numXis,
                         const int numEtas) {
  GraphPlan plan;
  plan.equations.resize(numEtas);

  const int numProducts = static_cast<int>(productDesign.n_rows);
  plan.products.resize(numProducts);
  for (int p = 0; p < numProducts; ++p) {
    plan.products[p].lastEta = -1;
    for (int v = 0; v < static_cast<int>(productDesign.n_cols); ++v) {
      const int exponent = static_cast<int>(productDesign(p, v));
      if (exponent <= 0) continue;
      plan.products[p].factors.push_back({v, exponent});
      if (v >= numXis) plan.products[p].lastEta =
        std::max(plan.products[p].lastEta, v - numXis);
    }
  }

  for (int equation = 0; equation < numEtas; ++equation) {
    for (int x = 0; x < numXis; ++x)
      if (gammaXi(equation, x) != 0.0)
        plan.equations[equation].xiLinear.push_back({x, -1, gammaXi(equation, x)});
    for (int eta = 0; eta < numEtas; ++eta)
      if (gammaEta(equation, eta) != 0.0)
        plan.equations[equation].etaLinear.push_back({eta, -1, gammaEta(equation, eta)});
    for (int p = 0; p < numProducts; ++p) {
      if (omega(equation, p) == 0.0) continue;
      // R validates evaluation order, but a product whose endogenous factors
      // are not yet available would silently read zeros here, so guard it.
      if (plan.products[p].lastEta >= equation)
        Rcpp::stop("A product term references an endogenous variable that is "
                   "not evaluated before the equation using it.");
      plan.equations[equation].products.push_back({p, omega(equation, p)});
    }
  }
  return plan;
}


// Evaluate one product across all nodes. Exponents are applied by repeated
// multiplication over the nonzero factors rather than by any pooled matrix
// operation: a 0/1 design matrix used as a BLAS selector is what produced the
// `0 * -Inf = NaN` bug in the ordinal kernel.
arma::vec evaluateProduct(const ProductPlan& product, const arma::mat& xi,
                          const arma::mat& eta, const int numXis,
                          const arma::uword rows) {
  arma::vec value(rows, arma::fill::ones);
  for (const ProductFactor& factor : product.factors) {
    const arma::vec base = factor.variable < numXis ?
      xi.col(factor.variable) : eta.col(factor.variable - numXis);
    for (int power = 0; power < factor.exponent; ++power) value %= base;
  }
  return value;
}


// Evaluate the latent state matrix from transformed innovations. Products are
// cached across equations but can only be evaluated once their endogenous
// factors exist, so evaluation interleaves with the structural order rather
// than happening upfront.
arma::mat evaluateGraphStates(const GraphPlan& plan,
                              const arma::mat& innovations,
                              const arma::vec& beta0, const arma::vec& alpha,
                              const int numXis, const int numEtas,
                              arma::mat* productOut = nullptr) {
  const arma::uword rows = innovations.n_rows;
  arma::mat xi(rows, numXis, arma::fill::zeros);
  if (numXis > 0) {
    xi = innovations.cols(0, numXis - 1);
    xi.each_row() += beta0.t();
  }

  arma::mat eta(rows, numEtas, arma::fill::zeros);
  const arma::uword numProducts = plan.products.size();
  arma::mat productValue(rows, numProducts, arma::fill::zeros);
  std::vector<bool> evaluated(numProducts, false);

  for (int equation = 0; equation < numEtas; ++equation) {
    // Evaluate every product whose factors are now available, not only those
    // this equation references. A free coefficient currently sitting at zero is
    // absent from the plan, but its derivative still needs the product value,
    // so an unevaluated column would silently score as zero.
    for (arma::uword p = 0; p < numProducts; ++p)
      if (!evaluated[p] && plan.products[p].lastEta < equation) {
        productValue.col(p) = evaluateProduct(plan.products[p], xi, eta,
                                              numXis, rows);
        evaluated[p] = true;
      }

    arma::vec value = innovations.col(numXis + equation) + alpha(equation);
    for (const GraphTerm& term : plan.equations[equation].xiLinear)
      value += term.coefficient * xi.col(term.first);
    for (const GraphTerm& term : plan.equations[equation].etaLinear)
      value += term.coefficient * eta.col(term.first);
    for (const GraphProductTerm& term : plan.equations[equation].products)
      value += term.coefficient * productValue.col(term.product);
    eta.col(equation) = value;
  }

  // Every eta exists now, so anything still unevaluated can be finished.
  for (arma::uword p = 0; p < numProducts; ++p)
    if (!evaluated[p]) {
      productValue.col(p) = evaluateProduct(plan.products[p], xi, eta,
                                            numXis, rows);
      evaluated[p] = true;
    }
  if (productOut != nullptr) *productOut = productValue;
  return arma::join_rows(xi, eta);
}


// Derivative of one product with respect to one of its factors, across every
// node. Built by re-multiplying the remaining factors, never by dividing the
// product, so a zero-valued factor stays finite.
arma::vec productFactorDerivative(const ProductPlan& product,
                                  const std::size_t index,
                                  const arma::mat& xi, const arma::mat& eta,
                                  const int numXis, const arma::uword rows) {
  const ProductFactor& target = product.factors[index];
  const arma::vec targetBase = target.variable < numXis ?
    xi.col(target.variable) : eta.col(target.variable - numXis);

  arma::vec derivative(rows, arma::fill::ones);
  derivative *= static_cast<double>(target.exponent);
  for (int power = 0; power < target.exponent - 1; ++power)
    derivative %= targetBase;
  for (std::size_t j = 0; j < product.factors.size(); ++j) {
    if (j == index) continue;
    const ProductFactor& other = product.factors[j];
    const arma::vec base = other.variable < numXis ?
      xi.col(other.variable) : eta.col(other.variable - numXis);
    for (int power = 0; power < other.exponent; ++power) derivative %= base;
  }
  return derivative;
}


// A product at one state, together with its derivative with respect to each
// factor. `d/dx_v = e_v * x_v^(e_v - 1) * prod_{u != v} x_u^{e_u}` is built by
// re-multiplying the other factors rather than dividing the total, so a factor
// that happens to be zero cannot produce NaN. Products have very few factors,
// so the quadratic cost is irrelevant.
struct ProductDerivative {
  double value;
  std::vector<double> partial; // aligned with product.factors
};

ProductDerivative productAt(const ProductPlan& product,
                            const arma::rowvec& xi,
                            const arma::rowvec& eta,
                            const int numXis) {
  const std::size_t n = product.factors.size();
  ProductDerivative out;
  out.value = 1.0;
  out.partial.assign(n, 0.0);

  std::vector<double> base(n);
  for (std::size_t i = 0; i < n; ++i) {
    const ProductFactor& factor = product.factors[i];
    base[i] = factor.variable < numXis ? xi(factor.variable) :
      eta(factor.variable - numXis);
    out.value *= std::pow(base[i], factor.exponent);
  }

  for (std::size_t i = 0; i < n; ++i) {
    double derivative = product.factors[i].exponent *
      std::pow(base[i], product.factors[i].exponent - 1);
    for (std::size_t j = 0; j < n; ++j)
      if (j != i) derivative *= std::pow(base[j], product.factors[j].exponent);
    out.partial[i] = derivative;
  }
  return out;
}

double logLogisticCdf(const double x) {
  if (x >= 0.0) return -std::log1p(std::exp(-x));
  return x - std::log1p(std::exp(x));
}

// Category log-probabilities are consumed by 0/1-weighted sums -- a design
// matrix product in the common kernel, a dot product against posterior counts
// in the M-step. IEEE gives 0 * -Inf = NaN, so a single underflowed category
// at an extreme node would poison every observation at that node, not just the
// ones actually in that category. Representing an underflowed probability by
// the smallest usable log instead of -Inf keeps those zero-weighted terms at
// zero; exp() of it is still 0, so nothing else changes.
constexpr double LOG_PROBABILITY_FLOOR = -700.0;

double floorLogProbability(const double x) {
  return (std::isnan(x) || x < LOG_PROBABILITY_FLOOR) ?
    LOG_PROBABILITY_FLOOR : x;
}

double logDiffExp(const double a, const double b) {
  if (!std::isfinite(b) && b < 0.0) return a;
  if (b >= a) return R_NegInf;
  return a + std::log1p(-std::exp(b - a));
}

double ordinalLogProbability(const int code, const arma::rowvec& thresholds,
                             const double eta, const bool logistic) {
  const int categories = static_cast<int>(thresholds.n_elem) + 1;
  if (code < 1 || code > categories) return LOG_PROBABILITY_FLOOR;
  const double lower = code == 1 ? R_NegInf : thresholds(code - 2) - eta;
  const double upper = code == categories ? R_PosInf : thresholds(code - 1) - eta;
  if (!std::isfinite(upper) && upper > 0.0) {
    if (logistic) return floorLogProbability(logLogisticCdf(-lower));
    return floorLogProbability(R::pnorm(-lower, 0.0, 1.0, true, true));
  }
  const double logUpper = logistic ? logLogisticCdf(upper) :
    R::pnorm(upper, 0.0, 1.0, true, true);
  const double logLower = logistic ? logLogisticCdf(lower) :
    R::pnorm(lower, 0.0, 1.0, true, true);
  return floorLogProbability(logDiffExp(logUpper, logLower));
}

double logResponseDensity(const double x, const bool logistic) {
  if (logistic) return logLogisticCdf(x) + logLogisticCdf(-x);
  return R::dnorm(x, 0.0, 1.0, true);
}

struct OrdinalEvaluation {
  double logProbability;
  double lowerRatio;
  double upperRatio;
};

OrdinalEvaluation ordinalEvaluation(const int code,
                                    const arma::rowvec& thresholds,
                                    const double eta,
                                    const bool logistic) {
  const int categories = static_cast<int>(thresholds.n_elem) + 1;
  const double lower = code == 1 ? R_NegInf : thresholds(code - 2) - eta;
  const double upper = code == categories ? R_PosInf : thresholds(code - 1) - eta;
  if (logistic) {
    const double logLower = std::isfinite(lower) ?
      logLogisticCdf(lower) : R_NegInf;
    const double logUpper = std::isfinite(upper) ?
      logLogisticCdf(upper) : 0.0;
    const double logp = floorLogProbability(logDiffExp(logUpper, logLower));
    if (logp <= LOG_PROBABILITY_FLOOR) return {logp, 0.0, 0.0};
    return {
      logp,
      std::isfinite(lower) ?
        std::exp(logLower + logLogisticCdf(-lower) - logp) : 0.0,
      std::isfinite(upper) ?
        std::exp(logUpper + logLogisticCdf(-upper) - logp) : 0.0
    };
  }
  const double logp = ordinalLogProbability(code, thresholds, eta, false);
  if (logp <= LOG_PROBABILITY_FLOOR) return {logp, 0.0, 0.0};
  return {
    logp,
    std::isfinite(lower) ? std::exp(logResponseDensity(lower, false) - logp) : 0.0,
    std::isfinite(upper) ? std::exp(logResponseDensity(upper, false) - logp) : 0.0
  };
}

// The sign split is a branch, not a gather: the two `find()` calls plus the
// `elem()` scatter walked the vector four times and allocated two index
// vectors, where `logLogisticCdf()` already applies exactly these two formulas
// per element. Kept vectorised in name only -- callers are unchanged.
arma::vec logLogisticCdfVector(const arma::vec& x) {
  arma::vec out(x.n_elem);
  const double* in = x.memptr();
  double* destination = out.memptr();
  for (arma::uword i = 0; i < x.n_elem; ++i)
    destination[i] = logLogisticCdf(in[i]);
  return out;
}

struct OrdinalBlockEvaluation {
  arma::vec logProbability;
  arma::vec lowerRatio;
  arma::vec upperRatio;
};

OrdinalBlockEvaluation ordinalBlockEvaluation(const int code,
                                               const arma::rowvec& thresholds,
                                               const arma::vec& eta,
                                               const bool logistic) {
  const int categories = static_cast<int>(thresholds.n_elem) + 1;
  OrdinalBlockEvaluation out = {
    arma::vec(eta.n_elem), arma::vec(eta.n_elem, arma::fill::zeros),
    arma::vec(eta.n_elem, arma::fill::zeros)
  };
  if (!logistic) {
    for (arma::uword q = 0; q < eta.n_elem; ++q) {
      const OrdinalEvaluation value = ordinalEvaluation(
        code, thresholds, eta(q), false
      );
      out.logProbability(q) = value.logProbability;
      out.lowerRatio(q) = value.lowerRatio;
      out.upperRatio(q) = value.upperRatio;
    }
    return out;
  }
  arma::vec logLower, logUpper;
  if (code > 1) logLower = logLogisticCdfVector(thresholds(code - 2) - eta);
  if (code < categories) logUpper =
    logLogisticCdfVector(thresholds(code - 1) - eta);
  if (code == 1) {
    out.logProbability = logUpper;
  } else if (code == categories) {
    out.logProbability = logLogisticCdfVector(eta - thresholds(code - 2));
  } else {
    // Both bounds underflow to the same value once eta is extreme, and the
    // vectorised form would then evaluate (-Inf) - (-Inf) = NaN. Defer to the
    // scalar routine, which treats an underflowed category as probability zero.
    for (arma::uword q = 0; q < eta.n_elem; ++q)
      out.logProbability(q) = logDiffExp(logUpper(q), logLower(q));
  }
  out.logProbability.transform(floorLogProbability);

  // A category with zero probability has no well-defined threshold ratio;
  // dividing by it would yield Inf and poison the score. Zero is the value the
  // scalar path uses for a structurally absent bound.
  const arma::uvec degenerate =
    arma::find(out.logProbability <= LOG_PROBABILITY_FLOOR);
  if (code > 1) {
    const arma::vec lower = thresholds(code - 2) - eta;
    out.lowerRatio = arma::exp(
      logLower + logLogisticCdfVector(-lower) - out.logProbability
    );
    out.lowerRatio.elem(degenerate).zeros();
  }
  if (code < categories) {
    const arma::vec upper = thresholds(code - 1) - eta;
    out.upperRatio = arma::exp(
      logUpper + logLogisticCdfVector(-upper) - out.logProbability
    );
    out.upperRatio.elem(degenerate).zeros();
  }
  return out;
}

// Log-probabilities only, for the callers that never look at the threshold
// ratios: every likelihood kernel below reads `logProbability` and nothing
// else. Producing the ratios anyway costs two further vectorised CDF passes
// and two allocations per (observation, indicator) on the logit path -- about
// as much again as the probabilities themselves -- and the packed kernels ran
// that on 3.4 of the 5.3 full N-by-Q passes each EM iteration.
//
// The scalar `ordinalLogProbability()` agrees with the vectorised branches
// term for term: `logDiffExp(logUpper, -Inf)` returns `logUpper` for the first
// category, and the `upper == +Inf` branch is `logLogisticCdf(-lower)`, which
// is what the last category's `logLogisticCdfVector(eta - thresholds)` gives.
// Looping it therefore changes cost, not results.
//
// `destination` is strided so a row of a column-major N-by-Q kernel can be
// written in place, which also drops the `subvec` copy of the means.
void accumulateOrdinalLogProbability(const int code,
                                     const arma::rowvec& thresholds,
                                     const double* means,
                                     const arma::uword n,
                                     const bool logistic,
                                     double* destination,
                                     const arma::uword stride) {
  for (arma::uword q = 0; q < n; ++q)
    destination[q * stride] +=
      ordinalLogProbability(code, thresholds, means[q], logistic);
}

double weightedOrdinalLogProbability(const int code,
                                     const arma::rowvec& thresholds,
                                     const double* means,
                                     const double* weights,
                                     const arma::uword n,
                                     const bool logistic) {
  double total = 0.0;
  for (arma::uword q = 0; q < n; ++q)
    total += weights[q] *
      ordinalLogProbability(code, thresholds, means[q], logistic);
  return total;
}

struct ResponsePlan {
  std::vector<bool> ordered;
  std::vector<arma::uvec> thresholdUse;
  std::vector<arma::rowvec> thresholds;
};

ResponsePlan buildResponsePlan(const arma::uword indicators,
                               const arma::mat& thresholdMatrix,
                               const std::vector<int>& orderedIndex) {
  ResponsePlan plan;
  plan.ordered.assign(indicators, false);
  plan.thresholdUse.resize(indicators);
  plan.thresholds.resize(indicators);
  for (int index : orderedIndex)
    if (index >= 0 && static_cast<arma::uword>(index) < indicators)
      plan.ordered[index] = true;
  for (arma::uword column = 0; column < indicators; ++column) {
    if (!plan.ordered[column]) continue;
    plan.thresholdUse[column] = arma::find_finite(thresholdMatrix.row(column).t());
    plan.thresholds[column].set_size(plan.thresholdUse[column].n_elem);
    for (arma::uword k = 0; k < plan.thresholdUse[column].n_elem; ++k)
      plan.thresholds[column](k) = thresholdMatrix(column, plan.thresholdUse[column](k));
  }
  return plan;
}

ResponsePlan buildResponsePlan(const arma::uword indicators,
                               const arma::mat& thresholdMatrix,
                               const Rcpp::IntegerVector& orderedIndex) {
  return buildResponsePlan(indicators, thresholdMatrix,
                           std::vector<int>(orderedIndex.begin(),
                                            orderedIndex.end()));
}

struct GraphObservation {
  arma::rowvec values;
  std::vector<int> columns;
};

std::vector<GraphObservation> flattenGraphObservations(
  const Rcpp::List& dataR, const Rcpp::List& colidxR);

// The per-indicator category design matrices and the dense value matrix that
// used to live here existed only for the shared-rule kernel, which tabulated
// Q by categories log-probabilities once and expanded them to N by Q with a
// BLAS product. That path is gone: one rule per observation is now the only
// evaluation path, and a shared rule is read from the same Q rows via a zero
// stride rather than through a separate implementation.
struct GraphWorkspace {
  std::vector<GraphObservation> observations;
  std::vector<int> orderedIndex;
  std::vector<bool> orderedLookup;

  GraphWorkspace(const Rcpp::List& dataR, const Rcpp::List& colidxR,
                 const Rcpp::IntegerVector& ordered) :
    observations(flattenGraphObservations(dataR, colidxR)),
    orderedIndex(ordered.begin(), ordered.end()) {
      int indicators = 0;
      for (const GraphObservation& observation : observations)
        for (const int column : observation.columns)
          indicators = std::max(indicators, column + 1);
      orderedLookup.assign(indicators, false);
      for (const int indicator : orderedIndex)
        if (indicator >= 0 && indicator < indicators)
          orderedLookup[indicator] = true;
    }
};

struct PreparedGraph {
  int numXis;
  int numEtas;
  int dimension;
  arma::mat latentChol;
  arma::mat gammaXi;
  arma::mat gammaEta;
  arma::mat productDesign;
  arma::mat omega;
  arma::vec beta0;
  arma::vec alpha;
  arma::mat lambda;
  arma::vec tau;
  arma::mat theta;
  GraphPlan structuralPlan;
  std::vector<LoadingTerm> loadingPlan;
  ResponsePlan responsePlan;

  PreparedGraph(const Rcpp::List& matrices, const int xis, const int etas,
                const std::vector<int>& orderedIndex) :
    numXis(xis), numEtas(etas), dimension(xis + etas),
    gammaXi(Rcpp::as<arma::mat>(matrices["gammaXi"])),
    gammaEta(Rcpp::as<arma::mat>(matrices["gammaEta"])),
    productDesign(Rcpp::as<arma::mat>(matrices["productDesign"])),
    omega(Rcpp::as<arma::mat>(matrices["omega"])),
    beta0(Rcpp::as<arma::vec>(matrices["beta0"])),
    alpha(Rcpp::as<arma::vec>(matrices["alpha"])),
    lambda(Rcpp::as<arma::mat>(matrices["lambdaX"])),
    tau(Rcpp::as<arma::vec>(matrices["tauX"])),
    theta(Rcpp::as<arma::mat>(matrices["thetaDelta"])) {
    const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
    const arma::mat cross = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
    const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
    const arma::mat thresholds = Rcpp::as<arma::mat>(matrices["thresholds"]);
    arma::mat covariance(dimension, dimension, arma::fill::zeros);
    if (numXis > 0)
      covariance.submat(0, 0, numXis - 1, numXis - 1) = A * A.t();
    if (numEtas > 0)
      covariance.submat(numXis, numXis, dimension - 1, dimension - 1) = psi;
    if (numXis > 0 && numEtas > 0) {
      covariance.submat(numXis, 0, dimension - 1, numXis - 1) = cross;
      covariance.submat(0, numXis, numXis - 1, dimension - 1) = cross.t();
    }
    if (!arma::chol(latentChol, covariance))
      Rcpp::stop("The unified latent covariance matrix is not positive definite.");
    structuralPlan = buildGraphPlan(gammaXi, gammaEta, productDesign,
                                    omega, numXis, numEtas);
    loadingPlan = buildLoadingPlan(lambda);
    responsePlan = buildResponsePlan(lambda.n_rows, thresholds, orderedIndex);
  }

  arma::mat states(const arma::mat& nodes) const {
    return evaluateGraphStates(structuralPlan, nodes * latentChol, beta0,
                               alpha, numXis, numEtas);
  }

  double negativeLogPosterior(const arma::vec& node,
                              const GraphObservation& observation,
                              arma::vec* gradient = nullptr,
                              const bool logistic = true) const {
    const arma::mat nodeMatrix = node.t();
    const arma::mat state = states(nodeMatrix);
    const arma::mat means = measurementMeans(state, tau, loadingPlan);
    arma::rowvec meanBar(lambda.n_rows, arma::fill::zeros);
    double logLikelihood = 0.0;
    for (arma::uword j = 0; j < observation.values.n_elem; ++j) {
      const int column = observation.columns[j];
      if (responsePlan.ordered[column]) {
        const int code = static_cast<int>(observation.values(j));
        const OrdinalEvaluation evaluation = ordinalEvaluation(
          code, responsePlan.thresholds[column], means(0, column), logistic
        );
        logLikelihood += evaluation.logProbability;
        meanBar(column) += evaluation.lowerRatio - evaluation.upperRatio;
      } else {
        const double variance = theta(column, column);
        const double residual = observation.values(j) - means(0, column);
        logLikelihood += -0.5 * (std::log(2.0 * M_PI * variance) +
          residual * residual / variance);
        meanBar(column) += residual / variance;
      }
    }
    if (gradient != nullptr) {
      arma::rowvec stateBar(dimension, arma::fill::zeros);
      for (const LoadingTerm& term : loadingPlan)
        stateBar(term.latent) += term.coefficient * meanBar(term.indicator);
      arma::rowvec xiBar = numXis > 0 ? stateBar.cols(0, numXis - 1) :
        arma::rowvec();
      arma::rowvec etaBar = numEtas > 0 ? stateBar.cols(numXis, dimension - 1) :
        arma::rowvec();
      const arma::rowvec xi = numXis > 0 ? state.cols(0, numXis - 1) :
        arma::rowvec();
      const arma::rowvec eta = numEtas > 0 ? state.cols(numXis, dimension - 1) :
        arma::rowvec();
      arma::rowvec innovationBar(dimension, arma::fill::zeros);
      for (int equation = numEtas - 1; equation >= 0; --equation) {
        const double b = etaBar(equation);
        innovationBar(numXis + equation) += b;
        for (const GraphTerm& term : structuralPlan.equations[equation].etaLinear)
          etaBar(term.first) += term.coefficient * b;
        for (const GraphTerm& term : structuralPlan.equations[equation].xiLinear)
          xiBar(term.first) += term.coefficient * b;
        for (const GraphProductTerm& term :
               structuralPlan.equations[equation].products) {
          const ProductPlan& product = structuralPlan.products[term.product];
          const ProductDerivative derivative = productAt(product, xi, eta, numXis);
          for (std::size_t i = 0; i < product.factors.size(); ++i) {
            const int variable = product.factors[i].variable;
            const double contribution =
              term.coefficient * b * derivative.partial[i];
            if (variable < numXis) xiBar(variable) += contribution;
            else etaBar(variable - numXis) += contribution;
          }
        }
      }
      if (numXis > 0) innovationBar.cols(0, numXis - 1) += xiBar;
      const arma::vec logLikelihoodGradient =
        (innovationBar * latentChol.t()).t();
      *gradient = node - logLikelihoodGradient;
    }
    return -logLikelihood + 0.5 * arma::dot(node, node);
  }
};

std::vector<GraphObservation> flattenGraphObservations(
    const Rcpp::List& dataR, const Rcpp::List& colidxR) {
  std::vector<GraphObservation> observations;
  for (int pattern = 0; pattern < dataR.size(); ++pattern) {
    const arma::mat values = Rcpp::as<arma::mat>(dataR[pattern]);
    const Rcpp::IntegerVector columnsR = colidxR[pattern];
    std::vector<int> columns(columnsR.begin(), columnsR.end());
    for (arma::uword i = 0; i < values.n_rows; ++i)
      observations.push_back({values.row(i), columns});
  }
  return observations;
}

arma::uword graphObservationCount(const Rcpp::List& dataR) {
  arma::uword total = 0;
  for (int pattern = 0; pattern < dataR.size(); ++pattern)
    total += Rcpp::as<arma::mat>(dataR[pattern]).n_rows;
  return total;
}

arma::mat responseLogKernel(const arma::mat& means,
                            const arma::mat& theta,
                            const ResponsePlan& responsePlan,
                            const Rcpp::List& dataR,
                            const Rcpp::List& colidxR,
                            const arma::uword nodesPerObservation,
                            const bool logistic,
                            const int ncores) {
  const arma::uword observations = graphObservationCount(dataR);
  arma::mat out(observations, nodesPerObservation, arma::fill::zeros);
  arma::uword offset = 0;
  for (int pattern = 0; pattern < dataR.size(); ++pattern) {
    const arma::mat values = Rcpp::as<arma::mat>(dataR[pattern]);
    const Rcpp::IntegerVector columns = colidxR[pattern];
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(ncores) if(ncores > 1) schedule(static)
    #endif
    for (int iInt = 0; iInt < static_cast<int>(values.n_rows); ++iInt) {
      const arma::uword i = static_cast<arma::uword>(iInt);
      const arma::uword first = (offset + i) * nodesPerObservation;
      const arma::uword last = first + nodesPerObservation - 1;
      for (arma::uword j = 0; j < values.n_cols; ++j) {
        const int column = columns[j];
        if (responsePlan.ordered[column]) {
          const int code = static_cast<int>(values(i, j));
          accumulateOrdinalLogProbability(
            code, responsePlan.thresholds[column],
            means.colptr(column) + first, nodesPerObservation, logistic,
            out.memptr() + (offset + i), out.n_rows
          );
        } else {
          const double variance = theta(column, column);
          const arma::vec residual = values(i, j) -
            means.col(column).subvec(first, last);
          out.row(offset + i) += (-0.5 *
            (std::log(2.0 * M_PI * variance) +
             arma::square(residual) / variance)).t();
        }
      }
    }
    offset += values.n_rows;
  }
  return out;
}


// `meanStride` is how far apart consecutive observations' node blocks sit in
// `means`. It is Q for a rule per observation, and 0 for a rule shared by every
// observation -- the shared case then reads the same Q rows for all of them,
// so a common rule stays Q nodes instead of being replicated into an N * Q
// matrix. At N = 10000 and Q = 3375 that is the difference between 81 KB and
// 810 MB of nodes, which is why `adaptive = "none"` keeps working at large N.
double responseWeightedKernel(
    const arma::mat& means,
    const arma::mat& theta,
    const ResponsePlan& responsePlan,
    const std::vector<GraphObservation>& observations,
    const arma::mat& weights,
    const bool logistic,
    const int ncores,
    const arma::uword meanStride) {
  if (weights.n_rows != observations.size())
    Rcpp::stop("Complete-data weights have incompatible dimensions.");
  const arma::uword Q = weights.n_cols;
  const arma::uword required = meanStride == 0 ?
    Q : observations.size() * meanStride;
  if (means.n_rows != required)
    Rcpp::stop("Graph nodes and complete-data weights are incompatible.");
  double total = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores > 1) \
    schedule(static) reduction(+:total)
  #endif
  for (int iInt = 0; iInt < static_cast<int>(observations.size()); ++iInt) {
    const arma::uword i = static_cast<arma::uword>(iInt);
    const GraphObservation& observation = observations[i];
    const arma::uword first = i * meanStride, last = first + Q - 1;
    const arma::vec observationWeights = weights.row(i).t();
    for (arma::uword j = 0; j < observation.values.n_elem; ++j) {
      const int column = observation.columns[j];
      if (responsePlan.ordered[column]) {
        total += weightedOrdinalLogProbability(
          static_cast<int>(observation.values(j)),
          responsePlan.thresholds[column],
          means.colptr(column) + first, observationWeights.memptr(),
          Q, logistic
        );
      } else {
        const double variance = theta(column, column);
        const arma::vec residual = observation.values(j) -
          means.col(column).subvec(first, last);
        total += arma::dot(observationWeights, -0.5 *
          (std::log(2.0 * M_PI * variance) +
           arma::square(residual) / variance));
      }
    }
  }
  return total;
}

// See `responseWeightedKernel` for what `meanStride` means: Q per observation,
// or 0 when one rule is shared by all of them.
arma::mat responseLogKernel(const arma::mat& means,
                            const arma::mat& theta,
                            const ResponsePlan& responsePlan,
                            const std::vector<GraphObservation>& observations,
                            const arma::uword nodesPerObservation,
                            const bool logistic,
                            const int ncores,
                            const arma::uword meanStride) {
  arma::mat out(observations.size(), nodesPerObservation, arma::fill::zeros);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores > 1) schedule(static)
  #endif
  for (int iInt = 0; iInt < static_cast<int>(observations.size()); ++iInt) {
    const arma::uword i = static_cast<arma::uword>(iInt);
    const GraphObservation& observation = observations[i];
    const arma::uword first = i * meanStride;
    const arma::uword last = first + nodesPerObservation - 1;
    for (arma::uword j = 0; j < observation.values.n_elem; ++j) {
      const int column = observation.columns[j];
      if (responsePlan.ordered[column]) {
        accumulateOrdinalLogProbability(
          static_cast<int>(observation.values(j)),
          responsePlan.thresholds[column],
          means.colptr(column) + first, nodesPerObservation, logistic,
          out.memptr() + i, out.n_rows
        );
      } else {
        const double variance = theta(column, column);
        const arma::vec residual = observation.values(j) -
          means.col(column).subvec(first, last);
        out.row(i) += (-0.5 * (std::log(2.0 * M_PI * variance) +
          arma::square(residual) / variance)).t();
      }
    }
  }
  return out;
}

arma::mat cholDirectional(const arma::mat& cholUpper, const arma::mat& dSigma) {
  const arma::mat left = arma::solve(arma::trimatl(cholUpper.t()), dSigma,
                                     arma::solve_opts::fast);
  arma::mat X = arma::solve(arma::trimatu(cholUpper), left.t(),
                            arma::solve_opts::fast).t();
  arma::mat phi(X.n_rows, X.n_cols, arma::fill::zeros);
  phi.elem(arma::trimatl_ind(arma::size(X), -1)) =
    X.elem(arma::trimatl_ind(arma::size(X), -1));
  phi.diag() = 0.5 * X.diag();
  return phi.t() * cholUpper;
}

// Quadrature weights are N * Q for a rule per observation and Q for a rule
// shared by every observation, matching how the nodes themselves are stored.
// The shared case is detected from the length rather than passed in: unlike
// the node matrix, the two shapes cannot collide, because Q == N * Q only when
// N == 1, where they are also the same thing.
Rcpp::List aggregateKernel(const arma::mat& logKernel,
                           const arma::vec& quadWeights,
                           const arma::vec& rowWeights) {
  const arma::uword observations = logKernel.n_rows;
  const arma::uword nodesPerObservation = logKernel.n_cols;
  const bool shared = quadWeights.n_elem == nodesPerObservation;
  if (!shared && quadWeights.n_elem != observations * nodesPerObservation)
    Rcpp::stop("Quadrature weights must have N * Q elements, or Q for a rule "
               "shared by every observation.");
  if (rowWeights.n_elem != observations)
    Rcpp::stop("Sampling weights must have one element per observation.");
  const arma::uword weightStride = shared ? 0 : nodesPerObservation;
  arma::mat posterior(logKernel.n_rows, logKernel.n_cols, arma::fill::zeros);
  arma::vec logDensity(logKernel.n_rows);
  double logLik = 0.0;
  for (arma::uword i = 0; i < logKernel.n_rows; ++i) {
    const arma::rowvec logWeights = arma::log(quadWeights.subvec(
      i * weightStride, i * weightStride + nodesPerObservation - 1
    )).t();
    const arma::rowvec joint = logKernel.row(i) + logWeights;
    const double maximum = joint.max();
    const double density = arma::accu(arma::exp(joint - maximum));
    logDensity(i) = maximum + std::log(density);
    posterior.row(i) = arma::exp(joint - logDensity(i));
    logLik += rowWeights(i) * logDensity(i);
  }
  return Rcpp::List::create(
    Rcpp::Named("logLik") = logLik,
    Rcpp::Named("posterior") = posterior,
    Rcpp::Named("logDensity") = logDensity
  );
}


arma::mat graphStates(const Rcpp::List& matrices, const arma::mat& nodes,
                      const int numXis, const int numEtas) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat cross = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat productDesign = Rcpp::as<arma::mat>(matrices["productDesign"]);
  const arma::mat omega = Rcpp::as<arma::mat>(matrices["omega"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);

  arma::mat latentCov(numXis + numEtas, numXis + numEtas, arma::fill::zeros);
  if (numXis > 0) latentCov.submat(0, 0, numXis - 1, numXis - 1) = A * A.t();
  if (numEtas > 0) {
    latentCov.submat(numXis, numXis, numXis + numEtas - 1,
                     numXis + numEtas - 1) = psi;
  }
  if (numXis > 0 && numEtas > 0) {
    latentCov.submat(numXis, 0, numXis + numEtas - 1, numXis - 1) = cross;
    latentCov.submat(0, numXis, numXis - 1, numXis + numEtas - 1) = cross.t();
  }

  arma::mat latentChol;
  if (!arma::chol(latentChol, latentCov))
    Rcpp::stop("The unified latent covariance matrix is not positive definite.");
  const GraphPlan plan = buildGraphPlan(
    gammaXi, gammaEta, productDesign, omega, numXis, numEtas
  );
  return evaluateGraphStates(plan, nodes * latentChol, beta0, alpha,
                             numXis, numEtas);
}

} // namespace


// [[Rcpp::export]]
SEXP lmsGraphWorkspaceCpp(const Rcpp::List& dataR,
                          const Rcpp::List& colidxR,
                          const Rcpp::IntegerVector& orderedIndex) {
  Rcpp::XPtr<GraphWorkspace> workspace(
    new GraphWorkspace(dataR, colidxR, orderedIndex), true
  );
  return workspace;
}










// Nodes per observation, given whether the rule is shared. `shared` is passed
// in rather than inferred from the row count: a shared rule whose Q happens to
// be a multiple of N (Q == N being the obvious case) is indistinguishable from
// a packed rule by shape alone, and guessing wrong silently reinterprets the
// whole quadrature. The caller always knows, from `rule$common`.
arma::uword graphNodesPerObservation(const arma::mat& nodes,
                                     const arma::uword observations,
                                     const bool shared) {
  if (observations == 0)
    Rcpp::stop("The graph workspace holds no observations.");
  if (shared) return nodes.n_rows;
  if (nodes.n_rows % observations != 0)
    Rcpp::stop("Packed graph nodes must have N * Q rows.");
  return nodes.n_rows / observations;
}


// [[Rcpp::export]]
arma::mat lmsGraphLogKernelWorkspaceCpp(const Rcpp::List& matrices,
                                        const arma::mat& nodes,
                                        const int numXis,
                                        const int numEtas,
                                        SEXP workspaceSEXP,
                                        const bool logistic = true,
                                        const int ncores = 1,
                                        const bool shared = false) {
  Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
  const PreparedGraph prepared(matrices, numXis, numEtas,
                               workspace->orderedIndex);
  const arma::mat states = prepared.states(nodes);
  const arma::mat means = measurementMeans(
    states, prepared.tau, prepared.loadingPlan
  );
  const arma::uword Q = graphNodesPerObservation(
    nodes, workspace->observations.size(), shared);
  return responseLogKernel(
    means, prepared.theta, prepared.responsePlan, workspace->observations,
    Q, logistic, ncores, shared ? 0 : Q
  );
}


// [[Rcpp::export]]
Rcpp::List lmsGraphPstepWorkspaceCpp(const Rcpp::List& matrices,
                                     const arma::mat& nodes,
                                     const int numXis,
                                     const int numEtas,
                                     SEXP workspaceSEXP,
                                     const arma::vec& quadWeights,
                                     const arma::vec& rowWeights,
                                     const bool logistic = true,
                                     const int ncores = 1,
                                     const bool shared = false) {
  const arma::mat logKernel = lmsGraphLogKernelWorkspaceCpp(
    matrices, nodes, numXis, numEtas, workspaceSEXP, logistic, ncores, shared
  );
  Rcpp::List out = aggregateKernel(logKernel, quadWeights, rowWeights);
  out["logKernel"] = logKernel;
  return out;
}


// [[Rcpp::export]]
double lmsGraphCompleteWorkspaceCpp(const Rcpp::List& matrices,
                                    const arma::mat& nodes,
                                    const int numXis,
                                    const int numEtas,
                                    SEXP workspaceSEXP,
                                    const arma::mat& weights,
                                    const bool logistic = true,
                                    const int ncores = 1,
                                    const bool shared = false) {
  Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
  const PreparedGraph prepared(matrices, numXis, numEtas,
                               workspace->orderedIndex);
  const arma::mat states = prepared.states(nodes);
  const arma::mat means = measurementMeans(
    states, prepared.tau, prepared.loadingPlan
  );
  const arma::uword Q = graphNodesPerObservation(
    nodes, workspace->observations.size(), shared);
  return responseWeightedKernel(
    means, prepared.theta, prepared.responsePlan, workspace->observations,
    weights, logistic, ncores, shared ? 0 : Q
  );
}


// [[Rcpp::export]]
arma::mat lmsGraphStatesCpp(const Rcpp::List& matrices,
                            const arma::mat& nodes,
                            const int numXis,
                            const int numEtas) {
  return graphStates(matrices, nodes, numXis, numEtas);
}


// [[Rcpp::export]]
Rcpp::List lmsGraphAdaptiveRuleCpp(const Rcpp::List& matrices,
                                   const arma::mat& baseNodes,
                                   const arma::vec& baseWeights,
                                   const arma::mat& starts,
                                   const int numXis,
                                   const int numEtas,
                                   const Rcpp::List& dataR,
                                   const Rcpp::List& colidxR,
                                   const Rcpp::IntegerVector& orderedIndex,
                                   const bool logistic = true,
                                   const int maxIterations = 25,
                                   const double tolerance = 1e-8,
                                   const double curvatureFloor = 1e-6,
                                   const double derivativeStep = 1e-4,
                                   const int ncores = 1,
                                   SEXP workspaceSEXP = R_NilValue) {
  std::vector<GraphObservation> temporaryObservations;
  const std::vector<GraphObservation>* observationData;
  std::vector<int> ordered;
  if (workspaceSEXP == R_NilValue) {
    temporaryObservations = flattenGraphObservations(dataR, colidxR);
    observationData = &temporaryObservations;
    ordered.assign(orderedIndex.begin(), orderedIndex.end());
  } else {
    Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
    observationData = &workspace->observations;
    ordered = workspace->orderedIndex;
  }
  const arma::uword observations = observationData->size();
  const arma::uword dimension = numXis + numEtas;
  const arma::uword Q = baseNodes.n_rows;
  if (baseNodes.n_cols != dimension || baseWeights.n_elem != Q)
    Rcpp::stop("The base quadrature rule has incompatible dimensions.");
  if (starts.n_rows != observations || starts.n_cols != dimension)
    Rcpp::stop("Adaptive quadrature starts must be an N by d matrix.");

  const PreparedGraph prepared(matrices, numXis, numEtas, ordered);
  arma::mat modes(observations, dimension, arma::fill::zeros);
  arma::cube transforms(dimension, dimension, observations, arma::fill::zeros);
  arma::mat packedNodes(observations * Q, dimension, arma::fill::zeros);
  arma::vec logWeights(observations * Q, arma::fill::zeros);
  arma::ivec convergence(observations, arma::fill::zeros);
  arma::uvec curvatureAdjusted(observations, arma::fill::zeros);
  const arma::vec logBasePrior = -0.5 * arma::sum(arma::square(baseNodes), 1) -
    0.5 * dimension * std::log(2.0 * M_PI);
  const arma::vec logBaseWeights = arma::log(baseWeights);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores > 1) schedule(dynamic)
  #endif
  for (int observationInt = 0;
       observationInt < static_cast<int>(observations); ++observationInt) {
      const arma::uword observation = static_cast<arma::uword>(observationInt);
      const GraphObservation& row = (*observationData)[observation];
      auto objective = [&](const arma::vec& z) {
        return prepared.negativeLogPosterior(z, row, nullptr, logistic);
      };
      auto objectiveGradient = [&](const arma::vec& z) {
        arma::vec value;
        prepared.negativeLogPosterior(z, row, &value, logistic);
        return value;
      };
      arma::vec mode = starts.row(observation).t();
      int code = 1;
      const double gradientTolerance = std::max(1e-6, std::sqrt(tolerance));
      arma::mat inverseHessian(dimension, dimension, arma::fill::eye);
      arma::vec gradient = objectiveGradient(mode);
      for (int iteration = 0; iteration < maxIterations; ++iteration) {
        if (!gradient.is_finite()) { code = 2; break; }
        if (arma::abs(gradient).max() < gradientTolerance) { code = 0; break; }
        arma::vec direction = -inverseHessian * gradient;
        double directionalDerivative = arma::dot(gradient, direction);
        if (!direction.is_finite() || directionalDerivative >= 0.0) {
          inverseHessian.eye();
          direction = -gradient;
          directionalDerivative = -arma::dot(gradient, gradient);
        }
        const double current = objective(mode);
        double scale = 1.0;
        bool accepted = false;
        arma::vec candidate;
        while (scale >= 1e-6) {
          candidate = mode + scale * direction;
          const double candidateValue = objective(candidate);
          if (std::isfinite(candidateValue) &&
              candidateValue <= current + 1e-4 * scale * directionalDerivative) {
            accepted = true; break;
          }
          scale *= 0.5;
        }
        if (!accepted) { code = 4; break; }
        const arma::vec newGradient = objectiveGradient(candidate);
        if (!newGradient.is_finite()) { code = 2; break; }
        const arma::vec stepTaken = candidate - mode;
        const arma::vec gradientChange = newGradient - gradient;
        const double curvature = arma::dot(gradientChange, stepTaken);
        if (curvature > 1e-10) {
          const double rho = 1.0 / curvature;
          const arma::mat identity(dimension, dimension, arma::fill::eye);
          const arma::mat left = identity - rho * stepTaken * gradientChange.t();
          inverseHessian = left * inverseHessian * left.t() +
            rho * stepTaken * stepTaken.t();
        } else {
          inverseHessian.eye();
        }
        mode = candidate;
        gradient = newGradient;
        if (arma::abs(stepTaken).max() < gradientTolerance) { code = 0; break; }
      }
      modes.row(observation) = mode.t();

      arma::vec eigenvalues;
      arma::mat hessian(dimension, dimension, arma::fill::zeros), eigenvectors;
      const arma::vec derivativeSteps = derivativeStep *
        (1.0 + arma::abs(mode));
      for (arma::uword j = 0; j < dimension; ++j) {
        arma::vec plus = mode, minus = mode;
        plus(j) += derivativeSteps(j); minus(j) -= derivativeSteps(j);
        hessian.col(j) = (objectiveGradient(plus) - objectiveGradient(minus)) /
          (2.0 * derivativeSteps(j));
      }
      hessian = 0.5 * (hessian + hessian.t());
      if (!arma::eig_sym(eigenvalues, eigenvectors, hessian) ||
          !eigenvalues.is_finite()) {
        eigenvalues.ones(dimension);
        eigenvectors.eye(dimension, dimension);
        curvatureAdjusted(observation) = 1;
      }
      for (arma::uword j = 0; j < dimension; ++j)
        if (eigenvalues(j) < curvatureFloor) {
          eigenvalues(j) = curvatureFloor;
          curvatureAdjusted(observation) = 1;
        }
      const arma::mat transform = eigenvectors *
        arma::diagmat(1.0 / arma::sqrt(eigenvalues));
      transforms.slice(observation) = transform;
      arma::mat adaptiveNodes = baseNodes * transform.t();
      adaptiveNodes.each_row() += mode.t();
      const arma::uword first = observation * Q, last = first + Q - 1;
      packedNodes.rows(first, last) = adaptiveNodes;
      const arma::vec logTargetPrior =
        -0.5 * arma::sum(arma::square(adaptiveNodes), 1) -
        0.5 * dimension * std::log(2.0 * M_PI);
      const double logJacobian = -0.5 * arma::accu(arma::log(eigenvalues));
      logWeights.subvec(first, last) = logBaseWeights + logTargetPrior +
        logJacobian - logBasePrior;
      convergence(observation) = code;
  }
  return Rcpp::List::create(
    Rcpp::Named("nodes") = packedNodes,
    Rcpp::Named("logWeights") = logWeights,
    Rcpp::Named("modes") = modes,
    Rcpp::Named("transforms") = transforms,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("curvatureAdjusted") = curvatureAdjusted
  );
}


// [[Rcpp::export]]
arma::mat lmsGraphLogKernelCpp(const Rcpp::List& matrices,
                               const arma::mat& nodes,
                               const int numXis,
                               const int numEtas,
                               const Rcpp::List& dataR,
                               const Rcpp::List& colidxR,
                               const Rcpp::IntegerVector& orderedIndex,
                               const bool logistic = true,
                               const int ncores = 1) {
  const arma::mat states = graphStates(matrices, nodes, numXis, numEtas);
  const arma::mat lambda = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::vec tau = Rcpp::as<arma::vec>(matrices["tauX"]);
  const arma::mat theta = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  const arma::mat thresholdMatrix = Rcpp::as<arma::mat>(matrices["thresholds"]);
  const std::vector<LoadingTerm> loadingPlan = buildLoadingPlan(lambda);
  const arma::mat means = measurementMeans(states, tau, loadingPlan);

  const ResponsePlan responsePlan = buildResponsePlan(
    lambda.n_rows, thresholdMatrix, orderedIndex
  );
  const arma::uword totalRows = graphObservationCount(dataR);
  if (totalRows == 0 || nodes.n_rows % totalRows != 0)
    Rcpp::stop("Packed graph nodes must have N * Q rows.");
  const arma::uword nodesPerObservation = nodes.n_rows / totalRows;
  return responseLogKernel(means, theta, responsePlan, dataR, colidxR,
                           nodesPerObservation, logistic, ncores);
}


// [[Rcpp::export]]
Rcpp::List lmsGraphAggregateCpp(const arma::mat& logKernel,
                                const arma::vec& quadWeights,
                                const arma::vec& rowWeights) {
  return aggregateKernel(logKernel, quadWeights, rowWeights);
}


// [[Rcpp::export]]
double lmsGraphWeightedKernelCpp(const arma::mat& logKernel,
                                 const arma::mat& weights) {
  if (logKernel.n_rows != weights.n_rows || logKernel.n_cols != weights.n_cols)
    Rcpp::stop("Likelihood kernel and posterior weights have incompatible dimensions.");
  return arma::accu(logKernel % weights);
}


// [[Rcpp::export]]
arma::vec lmsGraphScoreCpp(const Rcpp::List& matrices,
                           const arma::mat& nodes,
                           const int numXis,
                           const int numEtas,
                           const Rcpp::List& dataR,
                           const Rcpp::List& colidxR,
                           const Rcpp::IntegerVector& orderedIndex,
                           const arma::vec& quadWeights,
                           const arma::vec& rowWeights,
                           const arma::mat& completeWeights,
                           const bool observed,
                           const Rcpp::IntegerVector& block,
                           const Rcpp::IntegerVector& row,
                           const Rcpp::IntegerVector& col,
                           const Rcpp::LogicalVector& symmetric,
                           const bool logistic = true) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat cross = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat productDesign = Rcpp::as<arma::mat>(matrices["productDesign"]);
  const arma::mat omega = Rcpp::as<arma::mat>(matrices["omega"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);
  const arma::mat lambda = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::vec tau = Rcpp::as<arma::vec>(matrices["tauX"]);
  const arma::mat theta = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  const arma::mat thresholdMatrix = Rcpp::as<arma::mat>(matrices["thresholds"]);
  const arma::mat thresholdDelta = Rcpp::as<arma::mat>(matrices["thresholdDelta"]);
  const int latentDimension = numXis + numEtas;

  arma::mat latentCov(latentDimension, latentDimension, arma::fill::zeros);
  if (numXis > 0) latentCov.submat(0, 0, numXis - 1, numXis - 1) = A * A.t();
  if (numEtas > 0) latentCov.submat(numXis, numXis, latentDimension - 1,
                                    latentDimension - 1) = psi;
  if (numXis > 0 && numEtas > 0) {
    latentCov.submat(numXis, 0, latentDimension - 1, numXis - 1) = cross;
    latentCov.submat(0, numXis, numXis - 1, latentDimension - 1) = cross.t();
  }
  arma::mat latentChol;
  if (!arma::chol(latentChol, latentCov))
    Rcpp::stop("The unified latent covariance matrix is not positive definite.");
  const GraphPlan plan = buildGraphPlan(
    gammaXi, gammaEta, productDesign, omega, numXis, numEtas
  );
  arma::mat productValue;
  const arma::mat states = evaluateGraphStates(
    plan, nodes * latentChol, beta0, alpha, numXis, numEtas, &productValue
  );
  const arma::mat xi = numXis > 0 ? states.cols(0, numXis - 1) :
    arma::mat(nodes.n_rows, 0);
  const arma::mat eta = numEtas > 0 ?
    states.cols(numXis, numXis + numEtas - 1) : arma::mat(nodes.n_rows, 0);
  const std::vector<LoadingTerm> loadingPlan = buildLoadingPlan(lambda);
  const arma::mat means = measurementMeans(states, tau, loadingPlan);

  std::vector<bool> ordered(lambda.n_rows, false);
  for (int index : orderedIndex)
    if (index >= 0 && static_cast<arma::uword>(index) < lambda.n_rows)
      ordered[index] = true;

  const arma::mat logKernel = lmsGraphLogKernelCpp(
    matrices, nodes, numXis, numEtas, dataR, colidxR, orderedIndex, logistic
  );
  arma::mat likelihoodWeights;
  if (observed) {
    const Rcpp::List aggregate = aggregateKernel(logKernel, quadWeights, rowWeights);
    likelihoodWeights = Rcpp::as<arma::mat>(aggregate["posterior"]);
    likelihoodWeights.each_col() %= rowWeights;
  } else {
    likelihoodWeights = completeWeights;
  }
  const arma::uword nodesPerObservation = logKernel.n_cols;

  arma::vec score(block.size(), arma::fill::zeros);
  for (int parameter = 0; parameter < block.size(); ++parameter) {
    arma::mat dA(A.n_rows, A.n_cols, arma::fill::zeros);
    arma::mat dCross(cross.n_rows, cross.n_cols, arma::fill::zeros);
    arma::mat dPsi(psi.n_rows, psi.n_cols, arma::fill::zeros);
    arma::mat dGammaXi(gammaXi.n_rows, gammaXi.n_cols, arma::fill::zeros);
    arma::mat dGammaEta(gammaEta.n_rows, gammaEta.n_cols, arma::fill::zeros);
    arma::mat dOmega(omega.n_rows, omega.n_cols, arma::fill::zeros);
    arma::vec dBeta0(beta0.n_elem, arma::fill::zeros);
    arma::vec dAlpha(alpha.n_elem, arma::fill::zeros);
    arma::mat dLambda(lambda.n_rows, lambda.n_cols, arma::fill::zeros);
    arma::vec dTau(tau.n_elem, arma::fill::zeros);
    arma::mat dTheta(theta.n_rows, theta.n_cols, arma::fill::zeros);
    arma::mat dThreshold(thresholdMatrix.n_rows, thresholdMatrix.n_cols,
                         arma::fill::zeros);
    const int rr = row[parameter], cc = col[parameter], bb = block[parameter];
    switch (bb) {
    case 0: dLambda(rr, cc) = 1.0; break;
    case 2: dTau(rr) = 1.0; break;
    case 4:
      dTheta(rr, cc) = 1.0;
      if (symmetric[parameter] && rr != cc) dTheta(cc, rr) = 1.0;
      break;
    case 6: dA(rr, cc) = 1.0; break;
    case 7:
      dPsi(rr, cc) = 1.0;
      if (symmetric[parameter] && rr != cc) dPsi(cc, rr) = 1.0;
      break;
    case 8: dAlpha(rr) = 1.0; break;
    case 9: dBeta0(rr) = 1.0; break;
    case 10: dGammaXi(rr, cc) = 1.0; break;
    case 11: dGammaEta(rr, cc) = 1.0; break;
    case 19: dOmega(rr, cc) = 1.0; break;
    case 17: dCross(rr, cc) = 1.0; break;
    case 18: {
      arma::uvec use = arma::find_finite(thresholdDelta.row(rr).t());
      for (arma::uword position = 0; position < use.n_elem; ++position) {
        if (cc == 0 || static_cast<int>(use(position)) >= cc) {
          const double derivative = cc == 0 ? 1.0 :
            1.0 / (1.0 + std::exp(-thresholdDelta(rr, cc)));
          dThreshold(rr, use(position)) = derivative;
        }
      }
      break;
    }
    default: continue;
    }

    arma::mat dLatentCov(latentDimension, latentDimension, arma::fill::zeros);
    if (numXis > 0)
      dLatentCov.submat(0, 0, numXis - 1, numXis - 1) =
        dA * A.t() + A * dA.t();
    if (numEtas > 0)
      dLatentCov.submat(numXis, numXis, latentDimension - 1,
                        latentDimension - 1) = dPsi;
    if (numXis > 0 && numEtas > 0) {
      dLatentCov.submat(numXis, 0, latentDimension - 1, numXis - 1) = dCross;
      dLatentCov.submat(0, numXis, numXis - 1, latentDimension - 1) = dCross.t();
    }
    const arma::mat dInnovations = nodes * cholDirectional(latentChol, dLatentCov);
    arma::mat dXi(nodes.n_rows, numXis, arma::fill::zeros);
    if (numXis > 0) {
      dXi = dInnovations.cols(0, numXis - 1);
      dXi.each_row() += dBeta0.t();
    }
    arma::mat dEta(nodes.n_rows, numEtas, arma::fill::zeros);
    for (int j = 0; j < numEtas; ++j) {
      arma::vec dValue = dInnovations.col(numXis + j) + dAlpha(j);
      if (numXis > 0)
        dValue += dXi * gammaXi.row(j).t() + xi * dGammaXi.row(j).t();
      dValue += dEta * gammaEta.row(j).t() + eta * dGammaEta.row(j).t();
      // Forward mode over products: the coefficient's own direction plus the
      // chain rule through each factor. Every product is considered, not just
      // the ones in the plan, so a free coefficient currently at zero still
      // contributes its own direction.
      for (arma::uword pp = 0; pp < plan.products.size(); ++pp) {
        const ProductPlan& product = plan.products[pp];
        if (product.lastEta >= j) continue;
        if (omega(j, pp) != 0.0) {
          arma::vec dProduct(nodes.n_rows, arma::fill::zeros);
          for (std::size_t i = 0; i < product.factors.size(); ++i) {
            const int variable = product.factors[i].variable;
            const arma::vec dFactor = variable < numXis ? dXi.col(variable) :
              dEta.col(variable - numXis);
            dProduct += productFactorDerivative(product, i, xi, eta, numXis,
                                                nodes.n_rows) % dFactor;
          }
          dValue += omega(j, pp) * dProduct;
        }
        if (dOmega(j, pp) != 0.0)
          dValue += dOmega(j, pp) * productValue.col(pp);
      }
      dEta.col(j) = dValue;
    }
    const arma::mat dStates = arma::join_rows(dXi, dEta);
    arma::mat dMeans = dStates * lambda.t() + states * dLambda.t();
    dMeans.each_row() += dTau.t();

    arma::mat dLog(logKernel.n_rows, logKernel.n_cols, arma::fill::zeros);
    arma::uword offset = 0;
    for (int pattern = 0; pattern < dataR.size(); ++pattern) {
      const arma::mat values = Rcpp::as<arma::mat>(dataR[pattern]);
      const Rcpp::IntegerVector columns = colidxR[pattern];
      for (arma::uword i = 0; i < values.n_rows; ++i) {
        for (arma::uword q = 0; q < nodesPerObservation; ++q) {
          const arma::uword packed = (offset + i) * nodesPerObservation + q;
          double derivative = 0.0;
          for (arma::uword j = 0; j < values.n_cols; ++j) {
            const int column = columns[j];
            if (ordered[column]) {
              arma::uvec use = arma::find_finite(thresholdMatrix.row(column).t());
              const int code = static_cast<int>(values(i, j));
              const int categories = use.n_elem + 1;
              const double lower = code == 1 ? R_NegInf :
                thresholdMatrix(column, use(code - 2)) - means(packed, column);
              const double upper = code == categories ? R_PosInf :
                thresholdMatrix(column, use(code - 1)) - means(packed, column);
              arma::uvec itemUse = arma::find_finite(thresholdMatrix.row(column).t());
              arma::rowvec itemThresholds(itemUse.n_elem);
              for (arma::uword k = 0; k < itemUse.n_elem; ++k)
                itemThresholds(k) = thresholdMatrix(column, itemUse(k));
              const double logProbability = ordinalLogProbability(
                code, itemThresholds, means(packed, column), logistic
              );
              const double lowerRatio = std::isfinite(lower) ?
                std::exp(logResponseDensity(lower, logistic) - logProbability) : 0.0;
              const double upperRatio = std::isfinite(upper) ?
                std::exp(logResponseDensity(upper, logistic) - logProbability) : 0.0;
              derivative += (lowerRatio - upperRatio) * dMeans(packed, column);
              if (code > 1) derivative -= lowerRatio *
                dThreshold(column, use(code - 2));
              if (code < categories) derivative += upperRatio *
                dThreshold(column, use(code - 1));
            } else {
              const double variance = theta(column, column);
              const double residual = values(i, j) - means(packed, column);
              derivative += residual / variance * dMeans(packed, column);
              derivative += (-0.5 / variance + 0.5 * residual * residual /
                             (variance * variance)) * dTheta(column, column);
            }
          }
          dLog(offset + i, q) = derivative;
        }
      }
      offset += values.n_rows;
    }
    score(parameter) = arma::accu(likelihoodWeights % dLog);
  }
  return score;
}


// Reverse-mode score. The likelihood is evaluated once, response adjoints are
// accumulated over observations/nodes, and the structural graph is traversed
// once in reverse. Cholesky directions are evaluated only for covariance
// locations, never inside the observation/node loops.
// [[Rcpp::export]]
arma::vec lmsGraphReverseScoreCpp(const Rcpp::List& matrices,
                                  const arma::mat& nodes,
                                  const int numXis,
                                  const int numEtas,
                                  const Rcpp::List& dataR,
                                  const Rcpp::List& colidxR,
                                  const Rcpp::IntegerVector& orderedIndex,
                                  const arma::vec& quadWeights,
                                  const arma::vec& rowWeights,
                                  const arma::mat& completeWeights,
                                  const bool observed,
                                  const Rcpp::IntegerVector& block,
                                  const Rcpp::IntegerVector& row,
                                  const Rcpp::IntegerVector& col,
                                  const Rcpp::LogicalVector& symmetric,
                                  const bool logistic = true,
                                  const int ncores = 1,
                                  SEXP workspaceSEXP = R_NilValue,
                                  const bool shared = false) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat cross = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat productDesign = Rcpp::as<arma::mat>(matrices["productDesign"]);
  const arma::mat omega = Rcpp::as<arma::mat>(matrices["omega"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);
  const arma::mat lambda = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::vec tau = Rcpp::as<arma::vec>(matrices["tauX"]);
  const arma::mat theta = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  const arma::mat thresholdMatrix = Rcpp::as<arma::mat>(matrices["thresholds"]);
  const arma::mat thresholdDelta = Rcpp::as<arma::mat>(matrices["thresholdDelta"]);
  const int dimension = numXis + numEtas;

  arma::mat latentCov(dimension, dimension, arma::fill::zeros);
  if (numXis > 0) latentCov.submat(0, 0, numXis - 1, numXis - 1) = A * A.t();
  if (numEtas > 0)
    latentCov.submat(numXis, numXis, dimension - 1, dimension - 1) = psi;
  if (numXis > 0 && numEtas > 0) {
    latentCov.submat(numXis, 0, dimension - 1, numXis - 1) = cross;
    latentCov.submat(0, numXis, numXis - 1, dimension - 1) = cross.t();
  }
  arma::mat latentChol;
  if (!arma::chol(latentChol, latentCov))
    Rcpp::stop("The unified latent covariance matrix is not positive definite.");
  const GraphPlan plan = buildGraphPlan(
    gammaXi, gammaEta, productDesign, omega, numXis, numEtas
  );
  arma::mat productValue;
  const arma::mat states = evaluateGraphStates(
    plan, nodes * latentChol, beta0, alpha, numXis, numEtas, &productValue
  );
  const arma::mat xi = numXis > 0 ? states.cols(0, numXis - 1) :
    arma::mat(nodes.n_rows, 0);
  const arma::mat eta = numEtas > 0 ?
    states.cols(numXis, numXis + numEtas - 1) : arma::mat(nodes.n_rows, 0);
  const std::vector<LoadingTerm> loadingPlan = buildLoadingPlan(lambda);
  const arma::mat means = measurementMeans(states, tau, loadingPlan);
  const ResponsePlan responsePlan = buildResponsePlan(
    lambda.n_rows, thresholdMatrix, orderedIndex
  );
  std::vector<GraphObservation> temporaryObservations;
  const std::vector<GraphObservation>* observations;
  if (workspaceSEXP == R_NilValue) {
    temporaryObservations = flattenGraphObservations(dataR, colidxR);
    observations = &temporaryObservations;
  } else {
    Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
    observations = &workspace->observations;
  }
  const arma::uword totalRows = observations->size();
  const arma::uword nodesPerObservation =
    graphNodesPerObservation(nodes, totalRows, shared);
  const arma::uword meanStride = shared ? 0 : nodesPerObservation;
  arma::mat likelihoodWeights;
  if (observed) {
    const arma::mat logKernel = responseLogKernel(
      means, theta, responsePlan, *observations, nodesPerObservation,
      logistic, ncores, meanStride);
    const Rcpp::List aggregate = aggregateKernel(logKernel, quadWeights, rowWeights);
    likelihoodWeights = Rcpp::as<arma::mat>(aggregate["posterior"]);
    likelihoodWeights.each_col() %= rowWeights;
  } else {
    if (completeWeights.n_rows != totalRows ||
        completeWeights.n_cols != nodesPerObservation)
      Rcpp::stop("Complete-data weights have incompatible dimensions.");
    likelihoodWeights = completeWeights;
  }

  arma::mat meanBar(means.n_rows, means.n_cols, arma::fill::zeros);
  arma::mat thresholdBar(thresholdMatrix.n_rows, thresholdMatrix.n_cols,
                         arma::fill::zeros);
  arma::vec varianceBar(theta.n_rows, arma::fill::zeros);
  #ifdef _OPENMP
  #pragma omp parallel num_threads(ncores) if(ncores > 1)
  #endif
  {
    arma::mat thresholdLocal(thresholdBar.n_rows, thresholdBar.n_cols,
                             arma::fill::zeros);
    arma::vec varianceLocal(varianceBar.n_elem, arma::fill::zeros);
    // With a rule per observation each thread owns a disjoint block of
    // `meanBar` rows, so it can be written in place. With a shared rule every
    // observation writes the SAME Q rows, which would be a data race, so each
    // thread accumulates into its own copy and reduces below. That copy is
    // Q by J -- small, because a shared rule is exactly the case where the
    // node set is not replicated.
    const bool ownRows = !shared;
    arma::mat meanLocal;
    if (!ownRows) meanLocal.zeros(means.n_rows, means.n_cols);
    arma::mat& meanTarget = ownRows ? meanBar : meanLocal;
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int iInt = 0; iInt < static_cast<int>(totalRows); ++iInt) {
      const arma::uword i = static_cast<arma::uword>(iInt);
      const GraphObservation& observation = (*observations)[i];
      const arma::uword first = i * meanStride;
      const arma::uword last = first + nodesPerObservation - 1;
      const arma::vec weights = likelihoodWeights.row(i).t();
      for (arma::uword j = 0; j < observation.values.n_elem; ++j) {
          const int column = observation.columns[j];
          if (responsePlan.ordered[column]) {
            const arma::uvec& use = responsePlan.thresholdUse[column];
            const int code = static_cast<int>(observation.values(j));
            const int categories = use.n_elem + 1;
            const OrdinalBlockEvaluation evaluation = ordinalBlockEvaluation(
              code, responsePlan.thresholds[column],
              means.col(column).subvec(first, last), logistic
            );
            meanTarget.col(column).subvec(first, last) += weights %
              (evaluation.lowerRatio - evaluation.upperRatio);
            if (code > 1)
              thresholdLocal(column, use(code - 2)) -=
                arma::dot(weights, evaluation.lowerRatio);
            if (code < categories)
              thresholdLocal(column, use(code - 1)) +=
                arma::dot(weights, evaluation.upperRatio);
          } else {
            const double variance = theta(column, column);
            const arma::vec residual = observation.values(j) -
              means.col(column).subvec(first, last);
            meanTarget.col(column).subvec(first, last) +=
              weights % residual / variance;
            varianceLocal(column) += arma::dot(weights,
              -0.5 / variance + 0.5 * arma::square(residual) /
                (variance * variance));
        }
      }
    }
    #ifdef _OPENMP
    #pragma omp critical(lms_graph_adjoint_reduce)
    #endif
    {
      thresholdBar += thresholdLocal;
      varianceBar += varianceLocal;
      if (!ownRows) meanBar += meanLocal;
    }
  }

  const arma::vec tauBar = arma::sum(meanBar, 0).t();
  arma::mat stateBar(states.n_rows, states.n_cols, arma::fill::zeros);
  for (const LoadingTerm& term : loadingPlan)
    stateBar.col(term.latent) += term.coefficient * meanBar.col(term.indicator);
  arma::mat xiBar = numXis > 0 ? stateBar.cols(0, numXis - 1) :
    arma::mat(nodes.n_rows, 0, arma::fill::zeros);
  arma::mat etaBar = numEtas > 0 ? stateBar.cols(numXis, dimension - 1) :
    arma::mat(nodes.n_rows, 0, arma::fill::zeros);
  arma::mat equationBar(nodes.n_rows, numEtas, arma::fill::zeros);
  arma::vec alphaBar(alpha.n_elem, arma::fill::zeros);
  arma::vec beta0Bar(beta0.n_elem, arma::fill::zeros);
  arma::mat innovationBar(nodes.n_rows, dimension, arma::fill::zeros);

  for (int j = numEtas - 1; j >= 0; --j) {
    const arma::vec b = etaBar.col(j);
    equationBar.col(j) = b;
    innovationBar.col(numXis + j) += b;
    alphaBar(j) += arma::accu(b);
    for (const GraphTerm& term : plan.equations[j].etaLinear)
      etaBar.col(term.first) += term.coefficient * b;
    if (numXis > 0) {
      for (const GraphTerm& term : plan.equations[j].xiLinear)
        xiBar.col(term.first) += term.coefficient * b;
    }
    for (const GraphProductTerm& term : plan.equations[j].products) {
      const ProductPlan& product = plan.products[term.product];
      for (std::size_t i = 0; i < product.factors.size(); ++i) {
        const int variable = product.factors[i].variable;
        const arma::vec contribution = term.coefficient * b %
          productFactorDerivative(product, i, xi, eta, numXis, nodes.n_rows);
        if (variable < numXis) xiBar.col(variable) += contribution;
        else etaBar.col(variable - numXis) += contribution;
      }
    }
  }
  if (numXis > 0) {
    innovationBar.cols(0, numXis - 1) += xiBar;
    beta0Bar += arma::sum(xiBar, 0).t();
  }
  const arma::mat cholBar = nodes.t() * innovationBar;

  arma::vec score(block.size(), arma::fill::zeros);
  for (int parameter = 0; parameter < block.size(); ++parameter) {
    const int bb = block[parameter], rr = row[parameter], cc = col[parameter];
    switch (bb) {
    case 0:
      score(parameter) = arma::dot(meanBar.col(rr), states.col(cc));
      break;
    case 2: score(parameter) = tauBar(rr); break;
    case 4: score(parameter) = rr == cc ? varianceBar(rr) : 0.0; break;
    case 8: score(parameter) = alphaBar(rr); break;
    case 9: score(parameter) = beta0Bar(rr); break;
    case 10:
      score(parameter) = arma::dot(equationBar.col(rr), xi.col(cc));
      break;
    case 11:
      score(parameter) = arma::dot(equationBar.col(rr), eta.col(cc));
      break;
    case 19: {
      // d(eta_j)/d(omega_{j,p}) is simply the value of product p, which the
      // forward pass already cached.
      if (cc < 0 || static_cast<arma::uword>(cc) >= productValue.n_cols) break;
      score(parameter) = arma::dot(equationBar.col(rr), productValue.col(cc));
      break;
    }
    case 18: {
      if (rr < 0 || cc < 0 || static_cast<arma::uword>(rr) >= thresholdDelta.n_rows ||
          static_cast<arma::uword>(cc) >= thresholdDelta.n_cols ||
          !std::isfinite(thresholdDelta(rr, cc))) {
        score(parameter) = 0.0;
        break;
      }
      const arma::uvec use = arma::find_finite(thresholdDelta.row(rr).t());
      double value = 0.0;
      for (arma::uword position = 0; position < use.n_elem; ++position)
        if (cc == 0 || static_cast<int>(use(position)) >= cc)
          value += thresholdBar(rr, use(position));
      if (cc > 0) value *= 1.0 / (1.0 + std::exp(-thresholdDelta(rr, cc)));
      score(parameter) = value;
      break;
    }
    case 6: case 7: case 17: {
      arma::mat dSigma(dimension, dimension, arma::fill::zeros);
      if (bb == 6) {
        arma::mat dA(A.n_rows, A.n_cols, arma::fill::zeros); dA(rr, cc) = 1.0;
        dSigma.submat(0, 0, numXis - 1, numXis - 1) = dA * A.t() + A * dA.t();
      } else if (bb == 7) {
        arma::mat dPsi(psi.n_rows, psi.n_cols, arma::fill::zeros);
        dPsi(rr, cc) = 1.0;
        if (symmetric[parameter] && rr != cc) dPsi(cc, rr) = 1.0;
        dSigma.submat(numXis, numXis, dimension - 1, dimension - 1) = dPsi;
      } else {
        dSigma(numXis + rr, cc) = 1.0;
        dSigma(cc, numXis + rr) = 1.0;
      }
      score(parameter) = arma::accu(cholBar % cholDirectional(latentChol, dSigma));
      break;
    }
    default: break;
    }
  }
  return score;
}




// ---------------------------------------------------------------------------
// Structural sufficient statistics for the ECM M-step.
//
// PREMISE, NOT YET TRUE OF THE SHIPPED CODE (#28). Everything below assumes the
// E-step's nodes are held fixed through the M-step and read as latent values
// z = (xi, eta). They are not: `graphStates` treats them as standardised
// innovations, scaling by the latent Cholesky and pushing through
// beta0/alpha/Gamma/Omega, so z moves whenever a structural parameter does.
// That is why the ECM split is unwired -- see `lmsGraphBackend`. With z
// genuinely fixed the complete-data objective separates exactly:
//
//   Q(theta) = sum P log f(x | z; measurement)   <- needs the N by Q kernel
//            + sum P log p(z;   structural)      <- needs only the moments here
//
// and the structural half is a Gaussian log-density in
//   v = [xi - beta0 ; zeta],  zeta = eta - alpha - GammaXi xi - GammaEta eta
//                                    - Omega products(z),
// whose Jacobian with respect to eta is unit triangular in the recursive
// order, so it contributes nothing. Every term of v is AFFINE in the
// structural parameters once z is fixed, i.e. v = B(theta) h with
//   h = [1, xi, eta, products]
// fixed. Therefore
//   sum_iq P_iq v v' = B (sum_iq P_iq h h') B' = B S B'
// and one d-by-d matrix S, formed once per E-step, is all any structural
// evaluation needs. d = 1 + numXis + numEtas + numProducts, so the M-step
// structural objective is O(d^3) -- independent of both N and Q.
struct StructuralDesign {
  arma::uword dimension;   // rows of v
  arma::uword width;       // columns of h
};

StructuralDesign structuralDesign(const int numXis, const int numEtas,
                                  const arma::uword numProducts) {
  StructuralDesign out;
  out.dimension = static_cast<arma::uword>(numXis + numEtas);
  out.width = 1 + static_cast<arma::uword>(numXis + numEtas) + numProducts;
  return out;
}

// B maps h to v. Row i < numXis is xi_i - beta0_i; row numXis + m is
// eta_m minus its structural prediction.
arma::mat structuralMap(const Rcpp::List& matrices, const int numXis,
                        const int numEtas, const arma::uword numProducts) {
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat omega = Rcpp::as<arma::mat>(matrices["omega"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);
  const StructuralDesign design = structuralDesign(numXis, numEtas, numProducts);

  arma::mat B(design.dimension, design.width, arma::fill::zeros);
  for (int i = 0; i < numXis; ++i) {
    B(i, 0) = -beta0(i);
    B(i, 1 + i) = 1.0;
  }
  for (int m = 0; m < numEtas; ++m) {
    const arma::uword r = numXis + m;
    B(r, 0) = -alpha(m);
    for (int j = 0; j < numXis; ++j) B(r, 1 + j) = -gammaXi(m, j);
    for (int l = 0; l < numEtas; ++l)
      B(r, 1 + numXis + l) = (l == m ? 1.0 : 0.0) - gammaEta(m, l);
    for (arma::uword p = 0; p < numProducts; ++p)
      B(r, 1 + numXis + numEtas + p) = -omega(m, p);
  }
  return B;
}

arma::mat structuralCovariance(const Rcpp::List& matrices, const int numXis,
                               const int numEtas) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat cross = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::uword dimension = static_cast<arma::uword>(numXis + numEtas);
  arma::mat covariance(dimension, dimension, arma::fill::zeros);
  if (numXis > 0)
    covariance.submat(0, 0, numXis - 1, numXis - 1) = A * A.t();
  if (numEtas > 0)
    covariance.submat(numXis, numXis, dimension - 1, dimension - 1) = psi;
  if (numXis > 0 && numEtas > 0) {
    covariance.submat(numXis, 0, dimension - 1, numXis - 1) = cross;
    covariance.submat(0, numXis, numXis - 1, dimension - 1) = cross.t();
  }
  return covariance;
}


// S = sum_iq P_iq h h', with h built from the FIXED E-step nodes. `shared`
// has the same meaning as everywhere else: Q rows read by every observation,
// versus N * Q rows read one block each.
// [[Rcpp::export]]
Rcpp::List lmsGraphStructuralStatsCpp(const Rcpp::List& matrices,
                                      const arma::mat& nodes,
                                      const int numXis,
                                      const int numEtas,
                                      const arma::mat& posterior,
                                      const bool shared = false) {
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat productDesign = Rcpp::as<arma::mat>(matrices["productDesign"]);
  const arma::mat omega = Rcpp::as<arma::mat>(matrices["omega"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);

  arma::mat latentChol;
  if (!arma::chol(latentChol, structuralCovariance(matrices, numXis, numEtas)))
    Rcpp::stop("The unified latent covariance matrix is not positive definite.");
  const GraphPlan plan = buildGraphPlan(gammaXi, gammaEta, productDesign,
                                        omega, numXis, numEtas);
  arma::mat productValue;
  const arma::mat states = evaluateGraphStates(
    plan, nodes * latentChol, beta0, alpha, numXis, numEtas, &productValue);

  const arma::uword numProducts = productValue.n_cols;
  const StructuralDesign design = structuralDesign(numXis, numEtas, numProducts);
  const arma::uword N = posterior.n_rows;
  const arma::uword Q = posterior.n_cols;
  const arma::uword stride = shared ? 0 : Q;
  if (states.n_rows != (shared ? Q : N * Q))
    Rcpp::stop("Structural statistics require nodes matching the posterior.");

  // Collapsing the posterior over observations first turns the shared case
  // into Q outer products instead of N * Q -- the statistic is identical,
  // because h depends only on the node.
  arma::vec nodeMass;
  if (shared) nodeMass = arma::sum(posterior, 0).t();

  arma::mat S(design.width, design.width, arma::fill::zeros);
  double W = 0.0;
  arma::rowvec h(design.width);
  h(0) = 1.0;
  const arma::uword rows = shared ? Q : N * Q;
  for (arma::uword r = 0; r < rows; ++r) {
    double weight;
    if (shared) {
      weight = nodeMass(r);
    } else {
      weight = posterior(r / stride, r % stride);
    }
    if (weight == 0.0) continue;
    for (arma::uword c = 0; c < states.n_cols; ++c) h(1 + c) = states(r, c);
    for (arma::uword p = 0; p < numProducts; ++p)
      h(1 + states.n_cols + p) = productValue(r, p);
    S += weight * (h.t() * h);
    W += weight;
  }
  return Rcpp::List::create(
    Rcpp::Named("S") = S,
    Rcpp::Named("W") = W,
    Rcpp::Named("numProducts") = static_cast<int>(numProducts)
  );
}


// Structural half of the complete-data objective, from S and W alone.
// [[Rcpp::export]]
double lmsGraphStructuralCompleteCpp(const Rcpp::List& matrices,
                                     const arma::mat& S,
                                     const double W,
                                     const int numXis,
                                     const int numEtas,
                                     const int numProducts) {
  const arma::mat B = structuralMap(matrices, numXis, numEtas,
                                    static_cast<arma::uword>(numProducts));
  const arma::mat covariance = structuralCovariance(matrices, numXis, numEtas);
  arma::mat covarianceChol;
  if (!arma::chol(covarianceChol, covariance))
    Rcpp::stop("The unified latent covariance matrix is not positive definite.");
  const double logDeterminant =
    2.0 * arma::accu(arma::log(covarianceChol.diag()));
  const arma::mat quadratic = B * S * B.t();
  const arma::mat solved = arma::solve(arma::trimatu(covarianceChol),
    arma::solve(arma::trimatl(covarianceChol.t()), quadratic));
  const double dimension = static_cast<double>(numXis + numEtas);
  return -0.5 * (W * (dimension * std::log(2.0 * M_PI) + logDeterminant) +
                 arma::trace(solved));
}


// ---------------------------------------------------------------------------
// Measurement Newton step for the ECM M-step.
//
// With the latent nodes held fixed, the measurement half of the complete-data
// objective separates exactly across indicators: indicator l depends only on
// its own intercept, loadings, residual variance and thresholds, and on the
// nodes -- never on another indicator's parameters. Its Hessian is therefore
// block diagonal, and one traversal yields the objective together with every
// block's gradient and Hessian.
//
// Blocks are returned in the coordinates
//     [ intercept , loadings(active latents) , variance | thresholds ]
// with `variance` present for a continuous indicator and `thresholds` for an
// ordered one. R maps those onto free parameters, which for thresholds means
// chaining through the `thresholdDelta` parameterisation.
//
// This replaces the ~5 full N x Q passes nlminb spent line-searching the
// block with one. The second derivatives are the ones verified against central
// differences before any of this was written; with
//   a = lowerRatio, b = upperRatio, g = a - b, and f' the density derivative:
//     d2/dm2               = (f'(up) - f'(lo))/p - g^2
//     d2/dtau_c^2          = f'(up)/p - b^2
//     d2/dtau_{c-1}^2      = -f'(lo)/p - a^2
//     d2/dtau_c dtau_{c-1} = a*b
//     d2/dm dtau_c         = -f'(up)/p - g*b
//     d2/dm dtau_{c-1}     = f'(lo)/p + g*a
// f'/p is read off the ratios rather than recomputed: for the probit
// f'(x)/p = -x * f(x)/p, and for the logit f'(x)/p = (1 - 2F(x)) * f(x)/p.
// [[Rcpp::export]]
Rcpp::List lmsGraphMeasurementNewtonCpp(const Rcpp::List& matrices,
                                        const arma::mat& nodes,
                                        const int numXis,
                                        const int numEtas,
                                        const Rcpp::List& dataR,
                                        const Rcpp::List& colidxR,
                                        const Rcpp::IntegerVector& orderedIndex,
                                        const arma::mat& completeWeights,
                                        const Rcpp::List& activeLatents,
                                        const bool logistic = true,
                                        const int ncores = 1,
                                        SEXP workspaceSEXP = R_NilValue,
                                        const bool shared = false) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat cross = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat productDesign = Rcpp::as<arma::mat>(matrices["productDesign"]);
  const arma::mat omega = Rcpp::as<arma::mat>(matrices["omega"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);
  const arma::mat lambda = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::vec tau = Rcpp::as<arma::vec>(matrices["tauX"]);
  const arma::mat theta = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  const arma::mat thresholdMatrix = Rcpp::as<arma::mat>(matrices["thresholds"]);
  const int dimension = numXis + numEtas;

  arma::mat latentCov(dimension, dimension, arma::fill::zeros);
  if (numXis > 0) latentCov.submat(0, 0, numXis - 1, numXis - 1) = A * A.t();
  if (numEtas > 0)
    latentCov.submat(numXis, numXis, dimension - 1, dimension - 1) = psi;
  if (numXis > 0 && numEtas > 0) {
    latentCov.submat(numXis, 0, dimension - 1, numXis - 1) = cross;
    latentCov.submat(0, numXis, numXis - 1, dimension - 1) = cross.t();
  }
  arma::mat latentChol;
  if (!arma::chol(latentChol, latentCov))
    Rcpp::stop("The unified latent covariance matrix is not positive definite.");
  const GraphPlan plan = buildGraphPlan(
    gammaXi, gammaEta, productDesign, omega, numXis, numEtas
  );
  const arma::mat states = evaluateGraphStates(
    plan, nodes * latentChol, beta0, alpha, numXis, numEtas, nullptr
  );
  const std::vector<LoadingTerm> loadingPlan = buildLoadingPlan(lambda);
  const arma::mat means = measurementMeans(states, tau, loadingPlan);
  const ResponsePlan responsePlan = buildResponsePlan(
    lambda.n_rows, thresholdMatrix, orderedIndex
  );

  std::vector<GraphObservation> temporaryObservations;
  const std::vector<GraphObservation>* observations;
  if (workspaceSEXP == R_NilValue) {
    temporaryObservations = flattenGraphObservations(dataR, colidxR);
    observations = &temporaryObservations;
  } else {
    Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
    observations = &workspace->observations;
  }
  const arma::uword totalRows = observations->size();
  const arma::uword nodesPerObservation =
    graphNodesPerObservation(nodes, totalRows, shared);
  const arma::uword meanStride = shared ? 0 : nodesPerObservation;
  if (completeWeights.n_rows != totalRows ||
      completeWeights.n_cols != nodesPerObservation)
    Rcpp::stop("Complete-data weights have incompatible dimensions.");

  const arma::uword indicators = lambda.n_rows;

  // Per-indicator layout. `latents` are the latent columns whose loading on
  // this indicator can move; the intercept always occupies slot 0.
  std::vector<arma::uvec> latents(indicators);
  std::vector<arma::uword> meanSlots(indicators), totalSlots(indicators);
  for (arma::uword l = 0; l < indicators; ++l) {
    latents[l] = Rcpp::as<arma::uvec>(activeLatents[l]);
    meanSlots[l] = 1 + latents[l].n_elem;
    totalSlots[l] = meanSlots[l] +
      (responsePlan.ordered[l] ? responsePlan.thresholds[l].n_elem : 1);
  }

  std::vector<arma::vec> gradient(indicators);
  std::vector<arma::mat> hessian(indicators);
  for (arma::uword l = 0; l < indicators; ++l) {
    gradient[l].zeros(totalSlots[l]);
    hessian[l].zeros(totalSlots[l], totalSlots[l]);
  }
  double objective = 0.0;

  #ifdef _OPENMP
  #pragma omp parallel num_threads(ncores) if(ncores > 1)
  #endif
  {
    // The blocks are tiny -- (1 + loadings + categories)^2 per indicator -- so
    // every thread keeps its own copy and reduces once at the end. There is no
    // shared-rule race to worry about here for the same reason.
    std::vector<arma::vec> gradientLocal(indicators);
    std::vector<arma::mat> hessianLocal(indicators);
    for (arma::uword l = 0; l < indicators; ++l) {
      gradientLocal[l].zeros(totalSlots[l]);
      hessianLocal[l].zeros(totalSlots[l], totalSlots[l]);
    }
    double objectiveLocal = 0.0;
    arma::vec design;

    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int iInt = 0; iInt < static_cast<int>(totalRows); ++iInt) {
      const arma::uword i = static_cast<arma::uword>(iInt);
      const GraphObservation& observation = (*observations)[i];
      const arma::uword first = i * meanStride;
      const arma::uword last = first + nodesPerObservation - 1;
      const arma::vec weights = completeWeights.row(i).t();

      for (arma::uword j = 0; j < observation.values.n_elem; ++j) {
        const int column = observation.columns[j];
        const arma::uword slots = totalSlots[column];
        const arma::uword mSlots = meanSlots[column];
        const arma::uvec& use = latents[column];
        arma::vec& g = gradientLocal[column];
        arma::mat& H = hessianLocal[column];
        const arma::vec m = means.col(column).subvec(first, last);
        design.set_size(mSlots);
        design(0) = 1.0;

        if (responsePlan.ordered[column]) {
          const int code = static_cast<int>(observation.values(j));
          const arma::rowvec& thresholds = responsePlan.thresholds[column];
          const int categories = static_cast<int>(thresholds.n_elem) + 1;
          const OrdinalBlockEvaluation evaluation = ordinalBlockEvaluation(
            code, thresholds, m, logistic
          );
          objectiveLocal += arma::dot(weights, evaluation.logProbability);

          const arma::uword lowerSlot = mSlots + (code - 2);
          const arma::uword upperSlot = mSlots + (code - 1);

          for (arma::uword q = 0; q < weights.n_elem; ++q) {
            const double w = weights(q);
            if (w == 0.0) continue;
            const double a = evaluation.lowerRatio(q);
            const double b = evaluation.upperRatio(q);
            const double gm = a - b;

            // f'(x)/p from the ratio already in hand: no extra CDF evaluation
            // on the probit path, one logistic each on the logit path.
            double dLower = 0.0, dUpper = 0.0;
            if (code > 1) {
              const double x = thresholds(code - 2) - m(q);
              dLower = logistic ? a * (1.0 - 2.0 / (1.0 + std::exp(-x))) :
                                  -x * a;
            }
            if (code < categories) {
              const double x = thresholds(code - 1) - m(q);
              dUpper = logistic ? b * (1.0 - 2.0 / (1.0 + std::exp(-x))) :
                                  -x * b;
            }

            const double hmm = (dUpper - dLower) - gm * gm;
            for (arma::uword k = 1; k < mSlots; ++k)
              design(k) = states(first + q, use(k - 1));

            for (arma::uword r = 0; r < mSlots; ++r) {
              const double wd = w * design(r);
              g(r) += wd * gm;
              for (arma::uword c = 0; c <= r; ++c)
                H(r, c) += wd * hmm * design(c);
            }
            if (code > 1) {
              const double hll = -dLower - a * a;
              const double hml = dLower + gm * a;
              g(lowerSlot) -= w * a;
              H(lowerSlot, lowerSlot) += w * hll;
              for (arma::uword r = 0; r < mSlots; ++r)
                H(lowerSlot, r) += w * hml * design(r);
            }
            if (code < categories) {
              const double huu = dUpper - b * b;
              const double hmu = -dUpper - gm * b;
              g(upperSlot) += w * b;
              H(upperSlot, upperSlot) += w * huu;
              for (arma::uword r = 0; r < mSlots; ++r)
                H(upperSlot, r) += w * hmu * design(r);
            }
            if (code > 1 && code < categories)
              H(upperSlot, lowerSlot) += w * a * b;
          }

        } else {
          const double variance = theta(column, column);
          const arma::uword varSlot = slots - 1;
          const arma::vec residual = observation.values(j) - m;
          objectiveLocal += arma::dot(weights,
            -0.5 * std::log(2.0 * M_PI * variance) -
            0.5 * arma::square(residual) / variance);

          for (arma::uword q = 0; q < weights.n_elem; ++q) {
            const double w = weights(q);
            if (w == 0.0) continue;
            const double r0 = residual(q);
            const double gm = r0 / variance;
            const double hmm = -1.0 / variance;
            const double gv = -0.5 / variance +
              0.5 * r0 * r0 / (variance * variance);
            const double hvv = 0.5 / (variance * variance) -
              r0 * r0 / (variance * variance * variance);
            const double hmv = -r0 / (variance * variance);

            for (arma::uword k = 1; k < mSlots; ++k)
              design(k) = states(first + q, use(k - 1));

            for (arma::uword r = 0; r < mSlots; ++r) {
              const double wd = w * design(r);
              g(r) += wd * gm;
              for (arma::uword c = 0; c <= r; ++c)
                H(r, c) += wd * hmm * design(c);
            }
            g(varSlot) += w * gv;
            H(varSlot, varSlot) += w * hvv;
            for (arma::uword r = 0; r < mSlots; ++r)
              H(varSlot, r) += w * hmv * design(r);
          }
        }
      }
    }

    #ifdef _OPENMP
    #pragma omp critical(lms_graph_newton_reduce)
    #endif
    {
      objective += objectiveLocal;
      for (arma::uword l = 0; l < indicators; ++l) {
        gradient[l] += gradientLocal[l];
        hessian[l] += hessianLocal[l];
      }
    }
  }

  // Only the lower triangle was filled, since every contribution is symmetric.
  Rcpp::List gradientOut(indicators), hessianOut(indicators);
  for (arma::uword l = 0; l < indicators; ++l) {
    arma::mat& H = hessian[l];
    H = arma::symmatl(H);
    gradientOut[l] = gradient[l];
    hessianOut[l] = H;
  }
  return Rcpp::List::create(
    Rcpp::Named("objective") = objective,
    Rcpp::Named("gradient") = gradientOut,
    Rcpp::Named("hessian") = hessianOut
  );
}
