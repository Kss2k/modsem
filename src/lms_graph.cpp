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

struct ResponsePlan {
  std::vector<bool> ordered;
  std::vector<arma::uvec> thresholdUse;
  std::vector<arma::rowvec> thresholds;
};

ResponsePlan buildResponsePlan(const arma::uword indicators,
                               const arma::mat& thresholdMatrix,
                               const Rcpp::IntegerVector& orderedIndex) {
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
          for (arma::uword q = 0; q < nodesPerObservation; ++q)
            out(offset + i, q) += ordinalLogProbability(
              code, responsePlan.thresholds[column], means(first + q, column),
              logistic
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

double rowNegativeLogPosterior(const Rcpp::List& matrices,
                               const arma::vec& node,
                               const int numXis,
                               const int numEtas,
                               const arma::rowvec& values,
                               const Rcpp::IntegerVector& columns,
                               const ResponsePlan& responsePlan,
                               const bool logistic) {
  const arma::mat lambda = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::vec tau = Rcpp::as<arma::vec>(matrices["tauX"]);
  const arma::mat theta = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  const arma::mat nodeMatrix = node.t();
  const arma::mat states = graphStates(matrices, nodeMatrix, numXis, numEtas);
  const std::vector<LoadingTerm> loadingPlan = buildLoadingPlan(lambda);
  const arma::mat means = measurementMeans(states, tau, loadingPlan);
  double logLikelihood = 0.0;
  for (arma::uword j = 0; j < values.n_elem; ++j) {
    const int column = columns[j];
    if (responsePlan.ordered[column]) {
      logLikelihood += ordinalLogProbability(
        static_cast<int>(values(j)), responsePlan.thresholds[column],
        means(0, column), logistic
      );
    } else {
      const double variance = theta(column, column);
      const double residual = values(j) - means(0, column);
      logLikelihood += -0.5 * (std::log(2.0 * M_PI * variance) +
        residual * residual / variance);
    }
  }
  return -logLikelihood + 0.5 * arma::dot(node, node);
}

void numericalGradientHessian(const std::function<double(const arma::vec&)>& fn,
                              const arma::vec& x,
                              arma::vec& gradient,
                              arma::mat& hessian,
                              const double relativeStep) {
  const arma::uword dimension = x.n_elem;
  gradient.zeros(dimension);
  hessian.zeros(dimension, dimension);
  const double f0 = fn(x);
  arma::vec step = relativeStep * (1.0 + arma::abs(x));
  std::vector<double> plus(dimension), minus(dimension);
  for (arma::uword j = 0; j < dimension; ++j) {
    arma::vec xp = x, xm = x;
    xp(j) += step(j); xm(j) -= step(j);
    plus[j] = fn(xp); minus[j] = fn(xm);
    gradient(j) = (plus[j] - minus[j]) / (2.0 * step(j));
    hessian(j, j) = (plus[j] - 2.0 * f0 + minus[j]) /
      (step(j) * step(j));
  }
  for (arma::uword j = 0; j < dimension; ++j)
    for (arma::uword k = j + 1; k < dimension; ++k) {
      arma::vec pp = x, pm = x, mp = x, mm = x;
      pp(j) += step(j); pp(k) += step(k);
      pm(j) += step(j); pm(k) -= step(k);
      mp(j) -= step(j); mp(k) += step(k);
      mm(j) -= step(j); mm(k) -= step(k);
      const double value = (fn(pp) - fn(pm) - fn(mp) + fn(mm)) /
        (4.0 * step(j) * step(k));
      hessian(j, k) = hessian(k, j) = value;
    }
}

arma::vec numericalGradient(const std::function<double(const arma::vec&)>& fn,
                            const arma::vec& x,
                            const double relativeStep) {
  arma::vec gradient(x.n_elem, arma::fill::zeros);
  const arma::vec step = relativeStep * (1.0 + arma::abs(x));
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    arma::vec xp = x, xm = x;
    xp(j) += step(j); xm(j) -= step(j);
    gradient(j) = (fn(xp) - fn(xm)) / (2.0 * step(j));
  }
  return gradient;
}

} // namespace


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
                                   const double derivativeStep = 1e-4) {
  const arma::uword observations = graphObservationCount(dataR);
  const arma::uword dimension = numXis + numEtas;
  const arma::uword Q = baseNodes.n_rows;
  if (baseNodes.n_cols != dimension || baseWeights.n_elem != Q)
    Rcpp::stop("The base quadrature rule has incompatible dimensions.");
  if (starts.n_rows != observations || starts.n_cols != dimension)
    Rcpp::stop("Adaptive quadrature starts must be an N by d matrix.");

  const arma::mat lambda = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::mat thresholdMatrix = Rcpp::as<arma::mat>(matrices["thresholds"]);
  const ResponsePlan responsePlan = buildResponsePlan(
    lambda.n_rows, thresholdMatrix, orderedIndex
  );
  arma::mat modes(observations, dimension, arma::fill::zeros);
  arma::mat packedNodes(observations * Q, dimension, arma::fill::zeros);
  arma::vec logWeights(observations * Q, arma::fill::zeros);
  arma::ivec convergence(observations, arma::fill::zeros);
  arma::uvec curvatureAdjusted(observations, arma::fill::zeros);
  const arma::vec logBasePrior = -0.5 * arma::sum(arma::square(baseNodes), 1) -
    0.5 * dimension * std::log(2.0 * M_PI);
  const arma::vec logBaseWeights = arma::log(baseWeights);

  arma::uword observation = 0;
  for (int pattern = 0; pattern < dataR.size(); ++pattern) {
    const arma::mat values = Rcpp::as<arma::mat>(dataR[pattern]);
    const Rcpp::IntegerVector columns = colidxR[pattern];
    for (arma::uword i = 0; i < values.n_rows; ++i, ++observation) {
      if (observation % 100 == 0) Rcpp::checkUserInterrupt();
      const arma::rowvec rowValues = values.row(i);
      auto objective = [&](const arma::vec& z) {
        return rowNegativeLogPosterior(
          matrices, z, numXis, numEtas, rowValues, columns,
          responsePlan, logistic
        );
      };
      arma::vec mode = starts.row(observation).t();
      int code = 1;
      const double gradientTolerance = std::max(1e-6, std::sqrt(tolerance));
      arma::mat inverseHessian(dimension, dimension, arma::fill::eye);
      arma::vec gradient = numericalGradient(objective, mode, derivativeStep);
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
        const arma::vec newGradient = numericalGradient(
          objective, candidate, derivativeStep
        );
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

      arma::vec finalGradient, eigenvalues;
      arma::mat hessian, eigenvectors;
      numericalGradientHessian(objective, mode, finalGradient, hessian,
                               derivativeStep);
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
  }
  return Rcpp::List::create(
    Rcpp::Named("nodes") = packedNodes,
    Rcpp::Named("logWeights") = logWeights,
    Rcpp::Named("modes") = modes,
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
                                  const int ncores = 1) {
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
  const arma::uword totalRows = graphObservationCount(dataR);
  if (totalRows == 0 || nodes.n_rows % totalRows != 0)
    Rcpp::stop("Packed graph nodes must have N * Q rows.");
  const arma::uword nodesPerObservation = nodes.n_rows / totalRows;
  arma::mat likelihoodWeights;
  if (observed) {
    const arma::mat logKernel = responseLogKernel(
      means, theta, responsePlan, dataR, colidxR, nodesPerObservation,
      logistic, ncores
    );
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
  arma::uword offset = 0;
  for (int pattern = 0; pattern < dataR.size(); ++pattern) {
    const arma::mat values = Rcpp::as<arma::mat>(dataR[pattern]);
    const Rcpp::IntegerVector columns = colidxR[pattern];
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
      for (int iInt = 0; iInt < static_cast<int>(values.n_rows); ++iInt) {
        const arma::uword i = static_cast<arma::uword>(iInt);
        const arma::uword first = (offset + i) * nodesPerObservation;
        const arma::uword last = first + nodesPerObservation - 1;
        const arma::vec weights = likelihoodWeights.row(offset + i).t();
        for (arma::uword j = 0; j < values.n_cols; ++j) {
          const int column = columns[j];
          if (responsePlan.ordered[column]) {
            for (arma::uword q = 0; q < nodesPerObservation; ++q) {
              const arma::uword packed = first + q;
              const double weight = weights(q);
              const arma::uvec& use = responsePlan.thresholdUse[column];
              const int code = static_cast<int>(values(i, j));
              const int categories = use.n_elem + 1;
              const OrdinalEvaluation evaluation = ordinalEvaluation(
                code, responsePlan.thresholds[column], means(packed, column),
                logistic
              );
              meanBar(packed, column) += weight *
                (evaluation.lowerRatio - evaluation.upperRatio);
              if (code > 1)
                thresholdLocal(column, use(code - 2)) -=
                  weight * evaluation.lowerRatio;
              if (code < categories)
                thresholdLocal(column, use(code - 1)) +=
                  weight * evaluation.upperRatio;
            }
          } else {
            const double variance = theta(column, column);
            const arma::vec residual = values(i, j) -
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
    offset += values.n_rows;
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
