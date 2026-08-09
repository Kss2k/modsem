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

struct GraphEquationPlan {
  std::vector<GraphTerm> xiLinear;
  std::vector<GraphTerm> etaLinear;
  std::vector<GraphTerm> xiXi;
  std::vector<GraphTerm> xiEta;
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

std::vector<GraphEquationPlan> buildGraphPlan(const arma::mat& gammaXi,
                                               const arma::mat& gammaEta,
                                               const arma::mat& omegaXiXi,
                                               const arma::mat& omegaEtaXi,
                                               const int numXis,
                                               const int numEtas) {
  std::vector<GraphEquationPlan> plan(numEtas);
  for (int equation = 0; equation < numEtas; ++equation) {
    const int firstRow = equation * numXis;
    for (int x = 0; x < numXis; ++x) {
      if (gammaXi(equation, x) != 0.0)
        plan[equation].xiLinear.push_back({x, -1, gammaXi(equation, x)});
      for (int z = 0; z < numXis; ++z)
        if (omegaXiXi(firstRow + x, z) != 0.0)
          plan[equation].xiXi.push_back({x, z, omegaXiXi(firstRow + x, z)});
      for (int eta = 0; eta < numEtas; ++eta)
        if (omegaEtaXi(firstRow + x, eta) != 0.0)
          plan[equation].xiEta.push_back({x, eta, omegaEtaXi(firstRow + x, eta)});
    }
    for (int eta = 0; eta < numEtas; ++eta)
      if (gammaEta(equation, eta) != 0.0)
        plan[equation].etaLinear.push_back({eta, -1, gammaEta(equation, eta)});
  }
  return plan;
}

double logLogisticCdf(const double x) {
  if (x >= 0.0) return -std::log1p(std::exp(-x));
  return x - std::log1p(std::exp(x));
}

double logDiffExp(const double a, const double b) {
  if (!std::isfinite(b) && b < 0.0) return a;
  if (b >= a) return R_NegInf;
  return a + std::log1p(-std::exp(b - a));
}

double ordinalLogProbability(const int code, const arma::rowvec& thresholds,
                             const double eta, const bool logistic) {
  const int categories = static_cast<int>(thresholds.n_elem) + 1;
  if (code < 1 || code > categories) return R_NegInf;
  const double lower = code == 1 ? R_NegInf : thresholds(code - 2) - eta;
  const double upper = code == categories ? R_PosInf : thresholds(code - 1) - eta;
  if (!std::isfinite(upper) && upper > 0.0) {
    if (logistic) return logLogisticCdf(-lower);
    return R::pnorm(-lower, 0.0, 1.0, true, true);
  }
  const double logUpper = logistic ? logLogisticCdf(upper) :
    R::pnorm(upper, 0.0, 1.0, true, true);
  const double logLower = logistic ? logLogisticCdf(lower) :
    R::pnorm(lower, 0.0, 1.0, true, true);
  return logDiffExp(logUpper, logLower);
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
    const double logp = logDiffExp(logUpper, logLower);
    return {
      logp,
      std::isfinite(lower) ?
        std::exp(logLower + logLogisticCdf(-lower) - logp) : 0.0,
      std::isfinite(upper) ?
        std::exp(logUpper + logLogisticCdf(-upper) - logp) : 0.0
    };
  }
  const double logp = ordinalLogProbability(code, thresholds, eta, false);
  return {
    logp,
    std::isfinite(lower) ? std::exp(logResponseDensity(lower, false) - logp) : 0.0,
    std::isfinite(upper) ? std::exp(logResponseDensity(upper, false) - logp) : 0.0
  };
}

arma::vec logLogisticCdfVector(const arma::vec& x) {
  arma::vec out(x.n_elem);
  const arma::uvec nonnegative = arma::find(x >= 0.0);
  const arma::uvec negative = arma::find(x < 0.0);
  if (nonnegative.n_elem)
    out.elem(nonnegative) = -arma::log1p(arma::exp(-x.elem(nonnegative)));
  if (negative.n_elem)
    out.elem(negative) = x.elem(negative) -
      arma::log1p(arma::exp(x.elem(negative)));
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
    out.logProbability = logUpper +
      arma::log1p(-arma::exp(logLower - logUpper));
  }
  if (code > 1) {
    const arma::vec lower = thresholds(code - 2) - eta;
    out.lowerRatio = arma::exp(
      logLower + logLogisticCdfVector(-lower) - out.logProbability
    );
  }
  if (code < categories) {
    const arma::vec upper = thresholds(code - 1) - eta;
    out.upperRatio = arma::exp(
      logUpper + logLogisticCdfVector(-upper) - out.logProbability
    );
  }
  return out;
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

struct GraphIndicatorData {
  arma::uvec rows;
  arma::vec values;
  arma::mat categoryDesign;
};

std::vector<GraphObservation> flattenGraphObservations(
  const Rcpp::List& dataR, const Rcpp::List& colidxR);

struct GraphWorkspace {
  std::vector<GraphObservation> observations;
  std::vector<int> orderedIndex;
  std::vector<bool> orderedLookup;
  std::vector<GraphIndicatorData> indicators;
  arma::mat denseValues;

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

      std::vector<std::vector<arma::uword>> rowsByIndicator(indicators);
      std::vector<std::vector<double>> valuesByIndicator(indicators);
      for (arma::uword i = 0; i < observations.size(); ++i) {
        const GraphObservation& observation = observations[i];
        for (arma::uword j = 0; j < observation.values.n_elem; ++j) {
          const int indicator = observation.columns[j];
          rowsByIndicator[indicator].push_back(i);
          valuesByIndicator[indicator].push_back(observation.values(j));
        }
      }
      this->indicators.resize(indicators);
      denseValues.set_size(observations.size(), indicators);
      denseValues.fill(arma::datum::nan);
      for (int indicator = 0; indicator < indicators; ++indicator) {
        this->indicators[indicator].rows = arma::conv_to<arma::uvec>::from(
          rowsByIndicator[indicator]);
        this->indicators[indicator].values = arma::conv_to<arma::vec>::from(
          valuesByIndicator[indicator]);
        if (orderedLookup[indicator] && !valuesByIndicator[indicator].empty()) {
          const int categories = static_cast<int>(
            *std::max_element(valuesByIndicator[indicator].begin(),
                              valuesByIndicator[indicator].end()));
          this->indicators[indicator].categoryDesign.zeros(
            valuesByIndicator[indicator].size(), categories);
          for (arma::uword i = 0; i < valuesByIndicator[indicator].size(); ++i) {
            const int category = static_cast<int>(valuesByIndicator[indicator][i]) - 1;
            if (category >= 0 && category < categories)
              this->indicators[indicator].categoryDesign(i, category) = 1.0;
          }
        }
      }
      for (arma::uword i = 0; i < observations.size(); ++i)
        for (arma::uword j = 0; j < observations[i].values.n_elem; ++j)
          denseValues(i, observations[i].columns[j]) = observations[i].values(j);
    }
};

struct PreparedGraph {
  int numXis;
  int numEtas;
  int dimension;
  arma::mat latentChol;
  arma::mat gammaXi;
  arma::mat gammaEta;
  arma::mat omegaXiXi;
  arma::mat omegaEtaXi;
  arma::vec beta0;
  arma::vec alpha;
  arma::mat lambda;
  arma::vec tau;
  arma::mat theta;
  std::vector<GraphEquationPlan> structuralPlan;
  std::vector<LoadingTerm> loadingPlan;
  ResponsePlan responsePlan;

  PreparedGraph(const Rcpp::List& matrices, const int xis, const int etas,
                const std::vector<int>& orderedIndex) :
    numXis(xis), numEtas(etas), dimension(xis + etas),
    gammaXi(Rcpp::as<arma::mat>(matrices["gammaXi"])),
    gammaEta(Rcpp::as<arma::mat>(matrices["gammaEta"])),
    omegaXiXi(Rcpp::as<arma::mat>(matrices["omegaXiXi"])),
    omegaEtaXi(Rcpp::as<arma::mat>(matrices["omegaEtaXi"])),
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
    structuralPlan = buildGraphPlan(gammaXi, gammaEta, omegaXiXi,
                                    omegaEtaXi, numXis, numEtas);
    loadingPlan = buildLoadingPlan(lambda);
    responsePlan = buildResponsePlan(lambda.n_rows, thresholds, orderedIndex);
  }

  arma::mat states(const arma::mat& nodes) const {
    const arma::mat innovations = nodes * latentChol;
    arma::mat xi(nodes.n_rows, numXis, arma::fill::zeros);
    if (numXis > 0) {
      xi = innovations.cols(0, numXis - 1);
      xi.each_row() += beta0.t();
    }
    arma::mat eta(nodes.n_rows, numEtas, arma::fill::zeros);
    for (int equation = 0; equation < numEtas; ++equation) {
      arma::vec value = innovations.col(numXis + equation) + alpha(equation);
      for (const GraphTerm& term : structuralPlan[equation].xiLinear)
        value += term.coefficient * xi.col(term.first);
      for (const GraphTerm& term : structuralPlan[equation].etaLinear)
        value += term.coefficient * eta.col(term.first);
      for (const GraphTerm& term : structuralPlan[equation].xiXi)
        value += term.coefficient * xi.col(term.first) % xi.col(term.second);
      for (const GraphTerm& term : structuralPlan[equation].xiEta)
        value += term.coefficient * xi.col(term.first) % eta.col(term.second);
      eta.col(equation) = value;
    }
    return arma::join_rows(xi, eta);
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
        for (const GraphTerm& term : structuralPlan[equation].etaLinear)
          etaBar(term.first) += term.coefficient * b;
        for (const GraphTerm& term : structuralPlan[equation].xiLinear)
          xiBar(term.first) += term.coefficient * b;
        for (const GraphTerm& term : structuralPlan[equation].xiXi) {
          xiBar(term.first) += term.coefficient * b * xi(term.second);
          xiBar(term.second) += term.coefficient * b * xi(term.first);
        }
        for (const GraphTerm& term : structuralPlan[equation].xiEta) {
          xiBar(term.first) += term.coefficient * b * eta(term.second);
          etaBar(term.second) += term.coefficient * b * xi(term.first);
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
          const OrdinalBlockEvaluation evaluation = ordinalBlockEvaluation(
            code, responsePlan.thresholds[column],
            means.col(column).subvec(first, last), logistic
          );
          out.row(offset + i) += evaluation.logProbability.t();
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

// Response kernel for a quadrature rule shared by every observation.  The
// graph and measurement model are evaluated Q times, rather than N * Q times.
arma::mat responseLogKernelCommon(
    const arma::mat& means,
    const arma::mat& theta,
    const ResponsePlan& responsePlan,
    const GraphWorkspace& workspace,
    const bool logistic,
    const int ncores) {
  const std::vector<GraphObservation>& observations = workspace.observations;
  const arma::uword Q = means.n_rows;
  arma::mat out(observations.size(), Q, arma::fill::zeros);
  std::vector<arma::mat> ordinalTables(means.n_cols);
  for (arma::uword indicator = 0; indicator < means.n_cols; ++indicator) {
    const GraphIndicatorData& data = workspace.indicators[indicator];
    arma::mat contribution;
    if (responsePlan.ordered[indicator]) {
      const int categories = responsePlan.thresholds[indicator].n_elem + 1;
      ordinalTables[indicator].set_size(Q, categories);
      for (int category = 1; category <= categories; ++category)
        ordinalTables[indicator].col(category - 1) = ordinalBlockEvaluation(
          category, responsePlan.thresholds[indicator],
          means.col(indicator), logistic
        ).logProbability;
      contribution = data.categoryDesign * ordinalTables[indicator].t();
    } else {
      const double variance = theta(indicator, indicator);
      const arma::vec mu = means.col(indicator);
      contribution = data.values * mu.t() / variance;
      contribution.each_col() -= 0.5 * arma::square(data.values) / variance;
      contribution.each_row() -= 0.5 *
        (std::log(2.0 * M_PI * variance) + arma::square(mu) / variance).t();
    }
    if (data.rows.n_elem == observations.size() && data.rows(0) == 0 &&
        data.rows(data.rows.n_elem - 1) == observations.size() - 1) {
      out += contribution;
    } else {
      for (arma::uword i = 0; i < data.rows.n_elem; ++i)
        out.row(data.rows(i)) += contribution.row(i);
    }
  }
  return out;
}

double responseWeightedKernel(
    const arma::mat& means,
    const arma::mat& theta,
    const ResponsePlan& responsePlan,
    const std::vector<GraphObservation>& observations,
    const arma::mat& weights,
    const bool logistic,
    const int ncores) {
  if (weights.n_rows != observations.size())
    Rcpp::stop("Complete-data weights have incompatible dimensions.");
  const arma::uword Q = weights.n_cols;
  if (means.n_rows != observations.size() * Q)
    Rcpp::stop("Packed graph nodes and complete-data weights are incompatible.");
  double total = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores > 1) \
    schedule(static) reduction(+:total)
  #endif
  for (int iInt = 0; iInt < static_cast<int>(observations.size()); ++iInt) {
    const arma::uword i = static_cast<arma::uword>(iInt);
    const GraphObservation& observation = observations[i];
    const arma::uword first = i * Q, last = first + Q - 1;
    const arma::vec observationWeights = weights.row(i).t();
    for (arma::uword j = 0; j < observation.values.n_elem; ++j) {
      const int column = observation.columns[j];
      if (responsePlan.ordered[column]) {
        const OrdinalBlockEvaluation evaluation = ordinalBlockEvaluation(
          static_cast<int>(observation.values(j)),
          responsePlan.thresholds[column],
          means.col(column).subvec(first, last), logistic
        );
        total += arma::dot(observationWeights, evaluation.logProbability);
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

arma::mat responseLogKernel(const arma::mat& means,
                            const arma::mat& theta,
                            const ResponsePlan& responsePlan,
                            const std::vector<GraphObservation>& observations,
                            const arma::uword nodesPerObservation,
                            const bool logistic,
                            const int ncores) {
  arma::mat out(observations.size(), nodesPerObservation, arma::fill::zeros);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores > 1) schedule(static)
  #endif
  for (int iInt = 0; iInt < static_cast<int>(observations.size()); ++iInt) {
    const arma::uword i = static_cast<arma::uword>(iInt);
    const GraphObservation& observation = observations[i];
    const arma::uword first = i * nodesPerObservation;
    const arma::uword last = first + nodesPerObservation - 1;
    for (arma::uword j = 0; j < observation.values.n_elem; ++j) {
      const int column = observation.columns[j];
      if (responsePlan.ordered[column]) {
        const OrdinalBlockEvaluation evaluation = ordinalBlockEvaluation(
          static_cast<int>(observation.values(j)),
          responsePlan.thresholds[column],
          means.col(column).subvec(first, last), logistic
        );
        out.row(i) += evaluation.logProbability.t();
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

Rcpp::List aggregateKernel(const arma::mat& logKernel,
                           const arma::vec& quadWeights,
                           const arma::vec& rowWeights) {
  const arma::uword observations = logKernel.n_rows;
  const arma::uword nodesPerObservation = logKernel.n_cols;
  if (quadWeights.n_elem != observations * nodesPerObservation)
    Rcpp::stop("Packed quadrature weights must have N * Q elements.");
  if (rowWeights.n_elem != observations)
    Rcpp::stop("Sampling weights must have one element per observation.");
  arma::mat posterior(logKernel.n_rows, logKernel.n_cols, arma::fill::zeros);
  arma::vec logDensity(logKernel.n_rows);
  double logLik = 0.0;
  for (arma::uword i = 0; i < logKernel.n_rows; ++i) {
    const arma::rowvec logWeights = arma::log(quadWeights.subvec(
      i * nodesPerObservation, (i + 1) * nodesPerObservation - 1
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

Rcpp::List aggregateCommonKernel(const arma::mat& logKernel,
                                 const arma::vec& quadWeights,
                                 const arma::vec& rowWeights) {
  if (quadWeights.n_elem != logKernel.n_cols)
    Rcpp::stop("A common quadrature rule needs Q weights.");
  if (rowWeights.n_elem != logKernel.n_rows)
    Rcpp::stop("Sampling weights must have one element per observation.");
  const arma::rowvec logWeights = arma::log(quadWeights).t();
  arma::mat posterior(logKernel.n_rows, logKernel.n_cols, arma::fill::zeros);
  arma::vec logDensity(logKernel.n_rows);
  double logLik = 0.0;
  for (arma::uword i = 0; i < logKernel.n_rows; ++i) {
    const arma::rowvec joint = logKernel.row(i) + logWeights;
    const double maximum = joint.max();
    logDensity(i) = maximum + std::log(arma::accu(arma::exp(joint - maximum)));
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
  const arma::mat omegaXiXi = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat omegaEtaXi = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
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
  const arma::mat innovations = nodes * latentChol;
  arma::mat xi(nodes.n_rows, numXis, arma::fill::zeros);
  if (numXis > 0) {
    xi = innovations.cols(0, numXis - 1);
    xi.each_row() += beta0.t();
  }
  arma::mat eta(nodes.n_rows, numEtas, arma::fill::zeros);
  const std::vector<GraphEquationPlan> plan = buildGraphPlan(
    gammaXi, gammaEta, omegaXiXi, omegaEtaXi, numXis, numEtas
  );

  for (int j = 0; j < numEtas; ++j) {
    arma::vec value = innovations.col(numXis + j) + alpha(j);
    for (const GraphTerm& term : plan[j].xiLinear)
      value += term.coefficient * xi.col(term.first);
    for (const GraphTerm& term : plan[j].etaLinear)
      value += term.coefficient * eta.col(term.first);
    for (const GraphTerm& term : plan[j].xiXi)
      value += term.coefficient * xi.col(term.first) % xi.col(term.second);
    for (const GraphTerm& term : plan[j].xiEta)
      value += term.coefficient * xi.col(term.first) % eta.col(term.second);
    eta.col(j) = value;
  }
  return arma::join_rows(xi, eta);
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


// [[Rcpp::export]]
Rcpp::List lmsGraphCommonPstepCpp(const Rcpp::List& matrices,
                                  const arma::mat& nodes,
                                  const int numXis,
                                  const int numEtas,
                                  SEXP workspaceSEXP,
                                  const arma::vec& quadWeights,
                                  const arma::vec& rowWeights,
                                  const bool logistic = true,
                                  const int ncores = 1) {
  Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
  const PreparedGraph prepared(matrices, numXis, numEtas,
                               workspace->orderedIndex);
  const arma::mat states = prepared.states(nodes);
  const arma::mat means = measurementMeans(
    states, prepared.tau, prepared.loadingPlan
  );
  const arma::mat logKernel = responseLogKernelCommon(
    means, prepared.theta, prepared.responsePlan, *workspace,
    logistic, ncores
  );
  Rcpp::List out = aggregateCommonKernel(logKernel, quadWeights, rowWeights);
  out["logKernel"] = logKernel;
  return out;
}


// Fused common-rule E-step. Posterior probabilities exist only for one block
// at a time and are immediately reduced to M-step sufficient statistics.
// [[Rcpp::export]]
Rcpp::List lmsGraphCommonEstepCpp(const Rcpp::List& matrices,
                                  const arma::mat& nodes,
                                  const int numXis,
                                  const int numEtas,
                                  SEXP workspaceSEXP,
                                  const arma::vec& quadWeights,
                                  const arma::vec& rowWeights,
                                  const bool logistic = true,
                                  const int blockSize = 256,
                                  const int ncores = 1) {
  Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
  const PreparedGraph prepared(matrices, numXis, numEtas,
                               workspace->orderedIndex);
  const arma::mat states = prepared.states(nodes);
  const arma::mat means = measurementMeans(
    states, prepared.tau, prepared.loadingPlan);
  const arma::uword N = workspace->observations.size();
  const arma::uword Q = nodes.n_rows;
  const arma::uword J = means.n_cols;
  if (quadWeights.n_elem != Q || rowWeights.n_elem != N)
    Rcpp::stop("Common E-step quadrature or sampling weights are incompatible.");
  if (blockSize < 1) Rcpp::stop("E-step block size must be positive.");
  const arma::rowvec logQuadWeights = arma::log(quadWeights).t();

  int maxCategories = 1;
  std::vector<arma::mat> ordinalTables(J);
  for (arma::uword indicator = 0; indicator < J; ++indicator) {
    if (!prepared.responsePlan.ordered[indicator]) continue;
    const int categories = prepared.responsePlan.thresholds[indicator].n_elem + 1;
    maxCategories = std::max(maxCategories, categories);
    ordinalTables[indicator].set_size(categories, Q);
    for (int category = 1; category <= categories; ++category)
      ordinalTables[indicator].row(category - 1) = ordinalBlockEvaluation(
        category, prepared.responsePlan.thresholds[indicator],
        means.col(indicator), logistic).logProbability.t();
  }

  arma::mat mass(Q, J, arma::fill::zeros);
  arma::mat sum(Q, J, arma::fill::zeros);
  arma::mat sumsq(Q, J, arma::fill::zeros);
  arma::cube counts(Q, maxCategories, J, arma::fill::zeros);
  arma::vec logDensity(N);
  arma::vec nodeMass(Q, arma::fill::zeros);
  double logLik = 0.0;

  for (arma::uword first = 0; first < N; first += blockSize) {
    const arma::uword last = std::min(N - 1, first + blockSize - 1);
    const arma::uword B = last - first + 1;
    arma::mat posterior(B, Q, arma::fill::zeros);

    for (arma::uword indicator = 0; indicator < J; ++indicator) {
      const arma::vec values = workspace->denseValues.col(indicator).subvec(first, last);
      const arma::uvec use = arma::find_finite(values);
      if (use.n_elem == 0) continue;
      const arma::vec observedValues = values.elem(use);
      if (prepared.responsePlan.ordered[indicator]) {
        const arma::uvec category = arma::conv_to<arma::uvec>::from(observedValues - 1.0);
        const arma::mat contribution = ordinalTables[indicator].rows(category);
        if (use.n_elem == B) posterior += contribution;
        else posterior.rows(use) += contribution;
      } else {
        const double variance = prepared.theta(indicator, indicator);
        const arma::vec mu = means.col(indicator);
        arma::mat contribution = observedValues * mu.t() / variance;
        contribution.each_col() -= 0.5 * arma::square(observedValues) / variance;
        contribution.each_row() -= 0.5 *
          (std::log(2.0 * M_PI * variance) + arma::square(mu) / variance).t();
        if (use.n_elem == B) posterior += contribution;
        else posterior.rows(use) += contribution;
      }
    }

    for (arma::uword i = 0; i < B; ++i) {
      posterior.row(i) += logQuadWeights;
      const double maximum = posterior.row(i).max();
      posterior.row(i) = arma::exp(posterior.row(i) - maximum);
      const double density = arma::accu(posterior.row(i));
      logDensity(first + i) = maximum + std::log(density);
      logLik += rowWeights(first + i) * logDensity(first + i);
      posterior.row(i) *= rowWeights(first + i) / density;
      nodeMass += posterior.row(i).t();
    }

    for (arma::uword indicator = 0; indicator < J; ++indicator) {
      const arma::vec values = workspace->denseValues.col(indicator).subvec(first, last);
      const arma::uvec use = arma::find_finite(values);
      if (use.n_elem == 0) continue;
      const arma::vec observedValues = values.elem(use);
      const arma::mat weights = use.n_elem == B ? posterior : posterior.rows(use);
      mass.col(indicator) += arma::sum(weights, 0).t();
      if (prepared.responsePlan.ordered[indicator]) {
        const int categories = ordinalTables[indicator].n_rows;
        arma::mat design(use.n_elem, categories, arma::fill::zeros);
        for (arma::uword i = 0; i < use.n_elem; ++i)
          design(i, static_cast<int>(observedValues(i)) - 1) = 1.0;
        counts.slice(indicator).cols(0, categories - 1) += weights.t() * design;
      } else {
        sum.col(indicator) += weights.t() * observedValues;
        sumsq.col(indicator) += weights.t() * arma::square(observedValues);
      }
    }
  }
  (void)ncores;
  return Rcpp::List::create(
    Rcpp::Named("logLik") = logLik,
    Rcpp::Named("logDensity") = logDensity,
    Rcpp::Named("nodeMass") = nodeMass,
    Rcpp::Named("statistics") = Rcpp::List::create(
      Rcpp::Named("mass") = mass,
      Rcpp::Named("sum") = sum,
      Rcpp::Named("sumsq") = sumsq,
      Rcpp::Named("counts") = counts
    )
  );
}


// Reduce the weighted posterior to everything required by the conditionally
// independent continuous/ordinal measurement likelihood in the M-step.
// [[Rcpp::export]]
Rcpp::List lmsGraphSufficientStatsCpp(const arma::mat& posterior,
                                      SEXP workspaceSEXP,
                                      const int indicators,
                                      const int maxCategories,
                                      const int ncores = 1) {
  Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
  if (posterior.n_rows != workspace->observations.size())
    Rcpp::stop("Posterior rows and observations are incompatible.");
  const arma::uword Q = posterior.n_cols;
  arma::mat mass(Q, indicators, arma::fill::zeros);
  arma::mat sum(Q, indicators, arma::fill::zeros);
  arma::mat sumsq(Q, indicators, arma::fill::zeros);
  arma::cube counts(Q, std::max(1, maxCategories), indicators,
                    arma::fill::zeros);

  if (static_cast<int>(workspace->indicators.size()) > indicators)
    Rcpp::stop("The workspace has more indicators than the model matrices.");

  // Work indicator-wise. The posterior submatrix is column-major after the
  // row gather, allowing BLAS matrix-vector products to replace the former
  // Q * N * J scalar traversal. Each indicator owns its output column/slice.
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores > 1) schedule(dynamic)
  #endif
  for (int indicator = 0; indicator < indicators; ++indicator) {
    if (indicator >= static_cast<int>(workspace->indicators.size())) continue;
    const GraphIndicatorData& data = workspace->indicators[indicator];
    if (data.rows.n_elem == 0) continue;
    const bool complete = data.rows.n_elem == posterior.n_rows &&
      data.rows(0) == 0 && data.rows(data.rows.n_elem - 1) == posterior.n_rows - 1;
    arma::mat weights;
    if (!complete) weights = posterior.rows(data.rows);
    if (complete) mass.col(indicator) = arma::sum(posterior, 0).t();
    else mass.col(indicator) = arma::sum(weights, 0).t();
    const bool ordered = indicator <
      static_cast<int>(workspace->orderedLookup.size()) &&
      workspace->orderedLookup[indicator];
    if (ordered) {
      const arma::mat categoryCounts = complete ?
        posterior.t() * data.categoryDesign :
        weights.t() * data.categoryDesign;
      counts.slice(indicator).cols(0, categoryCounts.n_cols - 1) =
        categoryCounts;
    } else {
      if (complete) {
        sum.col(indicator) = posterior.t() * data.values;
        sumsq.col(indicator) = posterior.t() * arma::square(data.values);
      } else {
        sum.col(indicator) = weights.t() * data.values;
        sumsq.col(indicator) = weights.t() * arma::square(data.values);
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("mass") = mass,
    Rcpp::Named("sum") = sum,
    Rcpp::Named("sumsq") = sumsq,
    Rcpp::Named("counts") = counts
  );
}


// [[Rcpp::export]]
double lmsGraphSufficientCompleteCpp(const Rcpp::List& matrices,
                                     const arma::mat& nodes,
                                     const int numXis,
                                     const int numEtas,
                                     SEXP workspaceSEXP,
                                     const Rcpp::List& statistics,
                                     const bool logistic = true,
                                     const int ncores = 1) {
  Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
  const PreparedGraph prepared(matrices, numXis, numEtas,
                               workspace->orderedIndex);
  const arma::mat states = prepared.states(nodes);
  const arma::mat means = measurementMeans(
    states, prepared.tau, prepared.loadingPlan
  );
  const arma::mat mass = Rcpp::as<arma::mat>(statistics["mass"]);
  const arma::mat sum = Rcpp::as<arma::mat>(statistics["sum"]);
  const arma::mat sumsq = Rcpp::as<arma::mat>(statistics["sumsq"]);
  const arma::cube counts = Rcpp::as<arma::cube>(statistics["counts"]);
  if (mass.n_rows != nodes.n_rows || mass.n_cols != means.n_cols)
    Rcpp::stop("Graph sufficient statistics have incompatible dimensions.");
  double total = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores > 1) \
    schedule(static) reduction(+:total)
  #endif
  for (int indicatorInt = 0;
       indicatorInt < static_cast<int>(means.n_cols); ++indicatorInt) {
    const arma::uword indicator = static_cast<arma::uword>(indicatorInt);
    if (prepared.responsePlan.ordered[indicator]) {
      const int categories = prepared.responsePlan.thresholds[indicator].n_elem + 1;
      for (int category = 1; category <= categories; ++category) {
        const OrdinalBlockEvaluation evaluation = ordinalBlockEvaluation(
          category, prepared.responsePlan.thresholds[indicator],
          means.col(indicator), logistic
        );
        total += arma::dot(counts.slice(indicator).col(category - 1),
                           evaluation.logProbability);
      }
    } else {
      const double variance = prepared.theta(indicator, indicator);
      const arma::vec mu = means.col(indicator);
      const arma::vec rss = sumsq.col(indicator) -
        2.0 * mu % sum.col(indicator) + arma::square(mu) % mass.col(indicator);
      total += arma::accu(-0.5 * (mass.col(indicator) *
        std::log(2.0 * M_PI * variance) + rss / variance));
    }
  }
  return total;
}


// [[Rcpp::export]]
arma::mat lmsGraphLogKernelWorkspaceCpp(const Rcpp::List& matrices,
                                        const arma::mat& nodes,
                                        const int numXis,
                                        const int numEtas,
                                        SEXP workspaceSEXP,
                                        const bool logistic = true,
                                        const int ncores = 1) {
  Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
  const PreparedGraph prepared(matrices, numXis, numEtas,
                               workspace->orderedIndex);
  const arma::mat states = prepared.states(nodes);
  const arma::mat means = measurementMeans(
    states, prepared.tau, prepared.loadingPlan
  );
  const arma::uword observations = workspace->observations.size();
  if (observations == 0 || nodes.n_rows % observations != 0)
    Rcpp::stop("Packed graph nodes must have N * Q rows.");
  return responseLogKernel(
    means, prepared.theta, prepared.responsePlan, workspace->observations,
    nodes.n_rows / observations, logistic, ncores
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
                                     const int ncores = 1) {
  const arma::mat logKernel = lmsGraphLogKernelWorkspaceCpp(
    matrices, nodes, numXis, numEtas, workspaceSEXP, logistic, ncores
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
                                    const int ncores = 1) {
  Rcpp::XPtr<GraphWorkspace> workspace(workspaceSEXP);
  const PreparedGraph prepared(matrices, numXis, numEtas,
                               workspace->orderedIndex);
  const arma::mat states = prepared.states(nodes);
  const arma::mat means = measurementMeans(
    states, prepared.tau, prepared.loadingPlan
  );
  return responseWeightedKernel(
    means, prepared.theta, prepared.responsePlan, workspace->observations,
    weights, logistic, ncores
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
Rcpp::List lmsGraphPstepCpp(const Rcpp::List& matrices,
                            const arma::mat& nodes,
                            const int numXis,
                            const int numEtas,
                            const Rcpp::List& dataR,
                            const Rcpp::List& colidxR,
                            const Rcpp::IntegerVector& orderedIndex,
                            const arma::vec& quadWeights,
                            const arma::vec& rowWeights,
                            const bool logistic = true,
                            const int ncores = 1) {
  const arma::mat logKernel = lmsGraphLogKernelCpp(
    matrices, nodes, numXis, numEtas, dataR, colidxR, orderedIndex, logistic,
    ncores
  );
  Rcpp::List out = aggregateKernel(logKernel, quadWeights, rowWeights);
  out["logKernel"] = logKernel;
  return out;
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
  const arma::mat omegaXiXi = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat omegaEtaXi = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
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
  const arma::mat innovations = nodes * latentChol;
  arma::mat xi(nodes.n_rows, numXis, arma::fill::zeros);
  if (numXis > 0) {
    xi = innovations.cols(0, numXis - 1);
    xi.each_row() += beta0.t();
  }
  arma::mat eta(nodes.n_rows, numEtas, arma::fill::zeros);
  const std::vector<GraphEquationPlan> plan = buildGraphPlan(
    gammaXi, gammaEta, omegaXiXi, omegaEtaXi, numXis, numEtas
  );
  for (int j = 0; j < numEtas; ++j) {
    arma::vec value = innovations.col(numXis + j) + alpha(j);
    for (const GraphTerm& term : plan[j].xiLinear)
      value += term.coefficient * xi.col(term.first);
    for (const GraphTerm& term : plan[j].etaLinear)
      value += term.coefficient * eta.col(term.first);
    for (const GraphTerm& term : plan[j].xiXi)
      value += term.coefficient * xi.col(term.first) % xi.col(term.second);
    for (const GraphTerm& term : plan[j].xiEta)
      value += term.coefficient * xi.col(term.first) % eta.col(term.second);
    eta.col(j) = value;
  }
  const arma::mat states = arma::join_rows(xi, eta);
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
    arma::mat dOxx(omegaXiXi.n_rows, omegaXiXi.n_cols, arma::fill::zeros);
    arma::mat dOex(omegaEtaXi.n_rows, omegaEtaXi.n_cols, arma::fill::zeros);
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
    case 12: dOxx(rr, cc) = 1.0; break;
    case 13: dOex(rr, cc) = 1.0; break;
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
      if (numXis > 0) {
        const arma::uword first = j * numXis, last = first + numXis - 1;
        const arma::mat oxx = omegaXiXi.rows(first, last);
        const arma::mat oex = omegaEtaXi.rows(first, last);
        dValue += dXi * gammaXi.row(j).t() + xi * dGammaXi.row(j).t();
        dValue += arma::sum((dXi * oxx.t() + xi * dOxx.rows(first, last).t()) % xi +
                            (xi * oxx.t()) % dXi, 1);
        dValue += arma::sum((dXi * oex + xi * dOex.rows(first, last)) % eta +
                            (xi * oex) % dEta, 1);
      }
      dValue += dEta * gammaEta.row(j).t() + eta * dGammaEta.row(j).t();
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
                                  SEXP workspaceSEXP = R_NilValue) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat cross = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat omegaXiXi = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat omegaEtaXi = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
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
  const arma::mat innovations = nodes * latentChol;
  arma::mat xi(nodes.n_rows, numXis, arma::fill::zeros);
  if (numXis > 0) {
    xi = innovations.cols(0, numXis - 1);
    xi.each_row() += beta0.t();
  }
  arma::mat eta(nodes.n_rows, numEtas, arma::fill::zeros);
  const std::vector<GraphEquationPlan> plan = buildGraphPlan(
    gammaXi, gammaEta, omegaXiXi, omegaEtaXi, numXis, numEtas
  );
  for (int j = 0; j < numEtas; ++j) {
    arma::vec value = innovations.col(numXis + j) + alpha(j);
    for (const GraphTerm& term : plan[j].xiLinear)
      value += term.coefficient * xi.col(term.first);
    for (const GraphTerm& term : plan[j].etaLinear)
      value += term.coefficient * eta.col(term.first);
    for (const GraphTerm& term : plan[j].xiXi)
      value += term.coefficient * xi.col(term.first) % xi.col(term.second);
    for (const GraphTerm& term : plan[j].xiEta)
      value += term.coefficient * xi.col(term.first) % eta.col(term.second);
    eta.col(j) = value;
  }
  const arma::mat states = arma::join_rows(xi, eta);
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
  if (totalRows == 0 || nodes.n_rows % totalRows != 0)
    Rcpp::stop("Packed graph nodes must have N * Q rows.");
  const arma::uword nodesPerObservation = nodes.n_rows / totalRows;
  arma::mat likelihoodWeights;
  if (observed) {
    const arma::mat logKernel = responseLogKernel(
      means, theta, responsePlan, *observations, nodesPerObservation,
      logistic, ncores);
    const Rcpp::List aggregate = aggregateKernel(logKernel, quadWeights, rowWeights);
    likelihoodWeights = Rcpp::as<arma::mat>(aggregate["posterior"]);
    likelihoodWeights.each_col() %= rowWeights;
  } else {
    if (completeWeights.n_rows != totalRows ||
        completeWeights.n_cols != nodesPerObservation)
      Rcpp::stop("Complete-data weights have incompatible packed dimensions.");
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
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int iInt = 0; iInt < static_cast<int>(totalRows); ++iInt) {
      const arma::uword i = static_cast<arma::uword>(iInt);
      const GraphObservation& observation = (*observations)[i];
      const arma::uword first = i * nodesPerObservation;
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
            meanBar.col(column).subvec(first, last) += weights %
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
            meanBar.col(column).subvec(first, last) +=
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
    for (const GraphTerm& term : plan[j].etaLinear)
      etaBar.col(term.first) += term.coefficient * b;
    if (numXis > 0) {
      for (const GraphTerm& term : plan[j].xiLinear)
        xiBar.col(term.first) += term.coefficient * b;
      for (const GraphTerm& term : plan[j].xiXi) {
        xiBar.col(term.first) += term.coefficient * b % xi.col(term.second);
        xiBar.col(term.second) += term.coefficient * b % xi.col(term.first);
      }
      for (const GraphTerm& term : plan[j].xiEta) {
        xiBar.col(term.first) += term.coefficient * b % eta.col(term.second);
        etaBar.col(term.second) += term.coefficient * b % xi.col(term.first);
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
    case 12: {
      if (numXis <= 0) break;
      const int equation = rr / numXis;
      const int firstXi = rr % numXis;
      score(parameter) = arma::dot(
        equationBar.col(equation), xi.col(firstXi) % xi.col(cc)
      );
      break;
    }
    case 13: {
      if (numXis <= 0) break;
      const int equation = rr / numXis;
      const int firstXi = rr % numXis;
      score(parameter) = arma::dot(
        equationBar.col(equation), xi.col(firstXi) % eta.col(cc)
      );
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


// Reverse-mode complete-data score for one quadrature rule shared by a group.
// All observation-level work has already been reduced in the E-step.
// [[Rcpp::export]]
arma::vec lmsGraphSufficientScoreCpp(
    const Rcpp::List& matrices, const arma::mat& nodes,
    const int numXis, const int numEtas, const Rcpp::List& statistics,
    const Rcpp::IntegerVector& orderedIndex,
    const Rcpp::IntegerVector& block, const Rcpp::IntegerVector& row,
    const Rcpp::IntegerVector& col, const Rcpp::LogicalVector& symmetric,
    const bool logistic = true) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat cross = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat omegaXiXi = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat omegaEtaXi = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
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
  const arma::mat innovations = nodes * latentChol;
  arma::mat xi(nodes.n_rows, numXis, arma::fill::zeros);
  if (numXis > 0) { xi = innovations.cols(0, numXis - 1); xi.each_row() += beta0.t(); }
  arma::mat eta(nodes.n_rows, numEtas, arma::fill::zeros);
  const std::vector<GraphEquationPlan> plan = buildGraphPlan(
    gammaXi, gammaEta, omegaXiXi, omegaEtaXi, numXis, numEtas);
  for (int j = 0; j < numEtas; ++j) {
    arma::vec value = innovations.col(numXis + j) + alpha(j);
    for (const GraphTerm& term : plan[j].xiLinear)
      value += term.coefficient * xi.col(term.first);
    for (const GraphTerm& term : plan[j].etaLinear)
      value += term.coefficient * eta.col(term.first);
    for (const GraphTerm& term : plan[j].xiXi)
      value += term.coefficient * xi.col(term.first) % xi.col(term.second);
    for (const GraphTerm& term : plan[j].xiEta)
      value += term.coefficient * xi.col(term.first) % eta.col(term.second);
    eta.col(j) = value;
  }
  const arma::mat states = arma::join_rows(xi, eta);
  const std::vector<LoadingTerm> loadingPlan = buildLoadingPlan(lambda);
  const arma::mat means = measurementMeans(states, tau, loadingPlan);
  const ResponsePlan responsePlan = buildResponsePlan(
    lambda.n_rows, thresholdMatrix, orderedIndex);
  const arma::mat mass = Rcpp::as<arma::mat>(statistics["mass"]);
  const arma::mat sum = Rcpp::as<arma::mat>(statistics["sum"]);
  const arma::mat sumsq = Rcpp::as<arma::mat>(statistics["sumsq"]);
  const arma::cube counts = Rcpp::as<arma::cube>(statistics["counts"]);

  arma::mat meanBar(means.n_rows, means.n_cols, arma::fill::zeros);
  arma::mat thresholdBar(thresholdMatrix.n_rows, thresholdMatrix.n_cols,
                         arma::fill::zeros);
  arma::vec varianceBar(theta.n_rows, arma::fill::zeros);
  for (arma::uword indicator = 0; indicator < means.n_cols; ++indicator) {
    if (responsePlan.ordered[indicator]) {
      const arma::uvec& use = responsePlan.thresholdUse[indicator];
      const int categories = use.n_elem + 1;
      for (int category = 1; category <= categories; ++category) {
        const arma::vec weights = counts.slice(indicator).col(category - 1);
        const OrdinalBlockEvaluation evaluation = ordinalBlockEvaluation(
          category, responsePlan.thresholds[indicator],
          means.col(indicator), logistic);
        meanBar.col(indicator) += weights %
          (evaluation.lowerRatio - evaluation.upperRatio);
        if (category > 1)
          thresholdBar(indicator, use(category - 2)) -=
            arma::dot(weights, evaluation.lowerRatio);
        if (category < categories)
          thresholdBar(indicator, use(category - 1)) +=
            arma::dot(weights, evaluation.upperRatio);
      }
    } else {
      const double variance = theta(indicator, indicator);
      const arma::vec mu = means.col(indicator);
      meanBar.col(indicator) =
        (sum.col(indicator) - mu % mass.col(indicator)) / variance;
      const arma::vec rss = sumsq.col(indicator) - 2.0 * mu %
        sum.col(indicator) + arma::square(mu) % mass.col(indicator);
      varianceBar(indicator) = arma::accu(
        -0.5 * mass.col(indicator) / variance +
        0.5 * rss / (variance * variance));
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
    for (const GraphTerm& term : plan[j].etaLinear)
      etaBar.col(term.first) += term.coefficient * b;
    for (const GraphTerm& term : plan[j].xiLinear)
      xiBar.col(term.first) += term.coefficient * b;
    for (const GraphTerm& term : plan[j].xiXi) {
      xiBar.col(term.first) += term.coefficient * b % xi.col(term.second);
      xiBar.col(term.second) += term.coefficient * b % xi.col(term.first);
    }
    for (const GraphTerm& term : plan[j].xiEta) {
      xiBar.col(term.first) += term.coefficient * b % eta.col(term.second);
      etaBar.col(term.second) += term.coefficient * b % xi.col(term.first);
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
    case 0: score(parameter) = arma::dot(meanBar.col(rr), states.col(cc)); break;
    case 2: score(parameter) = tauBar(rr); break;
    case 4: score(parameter) = rr == cc ? varianceBar(rr) : 0.0; break;
    case 8: score(parameter) = alphaBar(rr); break;
    case 9: score(parameter) = beta0Bar(rr); break;
    case 10: score(parameter) = arma::dot(equationBar.col(rr), xi.col(cc)); break;
    case 11: score(parameter) = arma::dot(equationBar.col(rr), eta.col(cc)); break;
    case 12: {
      if (numXis > 0) {
        const int equation = rr / numXis, firstXi = rr % numXis;
        score(parameter) = arma::dot(equationBar.col(equation),
          xi.col(firstXi) % xi.col(cc));
      }
      break;
    }
    case 13: {
      if (numXis > 0) {
        const int equation = rr / numXis, firstXi = rr % numXis;
        score(parameter) = arma::dot(equationBar.col(equation),
          xi.col(firstXi) % eta.col(cc));
      }
      break;
    }
    case 18: {
      if (rr < 0 || cc < 0 || static_cast<arma::uword>(rr) >= thresholdDelta.n_rows ||
          static_cast<arma::uword>(cc) >= thresholdDelta.n_cols ||
          !std::isfinite(thresholdDelta(rr, cc))) break;
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
        dSigma(numXis + rr, cc) = 1.0; dSigma(cc, numXis + rr) = 1.0;
      }
      score(parameter) = arma::accu(cholBar % cholDirectional(latentChol, dSigma));
      break;
    }
    default: break;
    }
  }
  return score;
}
