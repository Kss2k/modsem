#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

inline double normal_cdf(const double x) {
  return R::pnorm(x, 0.0, 1.0, true, false);
}

inline double normal_density(const double x) {
  if (!std::isfinite(x)) return 0.0;
  return R::dnorm(x, 0.0, 1.0, false);
}

} // namespace


// Expected node-conditional log likelihood for a frozen E-step posterior.
// [[Rcpp::export]]
double lmsCatFrozenQCpp(const arma::mat& posterior,
                        const arma::mat& logKernel,
                        const arma::vec& rowWeights) {
  if (posterior.n_rows != logKernel.n_rows ||
      posterior.n_cols != logKernel.n_cols ||
      rowWeights.n_elem != posterior.n_rows)
    Rcpp::stop("Frozen-posterior Q inputs have incompatible dimensions.");
  return arma::dot(rowWeights, arma::sum(posterior % logKernel, 1));
}


// Fixed-quadrature cumulative-probit kernel for compressed response patterns.
// Codes are one-based, with NA or zero denoting a missing response.  `mu` is
// Q by J and is shared by the response patterns.  The returned predictor score
// is therefore Q by J; threshold scores follow the order of `thresholds`.
// [[Rcpp::export]]
Rcpp::List ordinalLogLikCpp(const Rcpp::IntegerMatrix& codes,
                            const arma::mat& mu,
                            const Rcpp::List& thresholds,
                            const arma::vec& quadWeights,
                            const arma::vec& logAdjust,
                            const arma::vec& patternWeights,
                            const bool gradient = false,
                            const bool details = false,
                            const bool posteriorDetails = true) {
  const arma::uword n = codes.nrow();
  const arma::uword jn = codes.ncol();
  const arma::uword qn = mu.n_rows;

  if (mu.n_cols != jn) Rcpp::stop("`mu` and `codes` have incompatible columns.");
  if (thresholds.size() != static_cast<R_xlen_t>(jn))
    Rcpp::stop("One threshold vector is required per indicator.");
  if (quadWeights.n_elem != qn || logAdjust.n_elem != qn)
    Rcpp::stop("Quadrature vectors must match the rows of `mu`.");
  if (patternWeights.n_elem != n)
    Rcpp::stop("Pattern weights must match the rows of `codes`.");

  std::vector<arma::vec> tau(jn);
  arma::uvec offsets(jn + 1, arma::fill::zeros);
  for (arma::uword j = 0; j < jn; ++j) {
    tau[j] = Rcpp::as<arma::vec>(thresholds[j]);
    offsets[j + 1] = offsets[j] + tau[j].n_elem;
  }

  arma::mat predictorScore(qn, jn, arma::fill::zeros);
  arma::vec thresholdScore(offsets[jn], arma::fill::zeros);
  arma::mat dmu(qn, jn, arma::fill::zeros);
  arma::mat upperDerivative(qn, jn, arma::fill::zeros);
  arma::mat lowerDerivative(qn, jn, arma::fill::zeros);
  arma::ivec category(jn, arma::fill::zeros);
  arma::vec logKernel(qn, arma::fill::zeros);
  arma::vec signedMass(qn, arma::fill::zeros);
  arma::mat posteriorOut(n, qn, arma::fill::zeros);
  arma::mat logKernelOut(n, qn, arma::fill::zeros);
  double logLik = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    logKernel.zeros();
    if (gradient) {
      dmu.zeros();
      upperDerivative.zeros();
      lowerDerivative.zeros();
    }

    for (arma::uword j = 0; j < jn; ++j) {
      const int code = codes(i, j);
      category[j] = code;
      if (code == NA_INTEGER || code <= 0) continue;
      const int categories = static_cast<int>(tau[j].n_elem) + 1;
      if (code > categories) Rcpp::stop("Category code exceeds its threshold count.");

      const double lo = code == 1 ? R_NegInf : tau[j][code - 2];
      const double hi = code == categories ? R_PosInf : tau[j][code - 1];
      for (arma::uword q = 0; q < qn; ++q) {
        const double eta = mu(q, j);
        const double Flo = std::isinf(lo) ? 0.0 : normal_cdf(lo - eta);
        const double Fhi = std::isinf(hi) ? 1.0 : normal_cdf(hi - eta);
        const double probability = std::max(Fhi - Flo, DBL_MIN);
        logKernel[q] += std::log(probability);
        if (gradient) {
          const double flo = normal_density(lo - eta);
          const double fhi = normal_density(hi - eta);
          dmu(q, j) = (flo - fhi) / probability;
          lowerDerivative(q, j) = -flo / probability;
          upperDerivative(q, j) = fhi / probability;
        }
      }
    }

    double scale = R_NegInf;
    for (arma::uword q = 0; q < qn; ++q) {
      if (quadWeights[q] == 0.0) continue;
      scale = std::max(scale, logKernel[q] + logAdjust[q] +
                       std::log(std::abs(quadWeights[q])));
    }
    if (!std::isfinite(scale))
      return Rcpp::List::create(Rcpp::Named("logLik") = R_NegInf);

    double integralScaled = 0.0;
    for (arma::uword q = 0; q < qn; ++q) {
      signedMass[q] = quadWeights[q] *
        std::exp(logKernel[q] + logAdjust[q] - scale);
      integralScaled += signedMass[q];
    }
    if (!(integralScaled > 0.0) || !std::isfinite(integralScaled))
      return Rcpp::List::create(Rcpp::Named("logLik") = R_NegInf);

    const double frequency = patternWeights[i];
    logLik += frequency * (scale + std::log(integralScaled));
    if (gradient || details) {
      logKernelOut.row(i) = logKernel.t();
      posteriorOut.row(i) = (signedMass / integralScaled).t();
    }
    if (!gradient) continue;

    for (arma::uword q = 0; q < qn; ++q) {
      const double posterior = frequency * signedMass[q] / integralScaled;
      for (arma::uword j = 0; j < jn; ++j) {
        const int code = category[j];
        if (code == NA_INTEGER || code <= 0) continue;
        predictorScore(q, j) += posterior * dmu(q, j);
        const int categories = static_cast<int>(tau[j].n_elem) + 1;
        if (code > 1)
          thresholdScore[offsets[j] + code - 2] += posterior * lowerDerivative(q, j);
        if (code < categories)
          thresholdScore[offsets[j] + code - 1] += posterior * upperDerivative(q, j);
      }
    }
  }

  if (!gradient && !details)
    return Rcpp::List::create(Rcpp::Named("logLik") = logLik);
  if (!gradient)
    return posteriorDetails ? Rcpp::List::create(
        Rcpp::Named("logLik") = logLik,
        Rcpp::Named("posterior") = posteriorOut,
        Rcpp::Named("log.kernel") = logKernelOut
      ) : Rcpp::List::create(
        Rcpp::Named("logLik") = logLik,
        Rcpp::Named("log.kernel") = logKernelOut
      );
  return Rcpp::List::create(
    Rcpp::Named("logLik") = logLik,
    Rcpp::Named("predictor.score") = predictorScore,
    Rcpp::Named("threshold.score") = thresholdScore,
    Rcpp::Named("posterior") = posteriorOut,
    Rcpp::Named("log.kernel") = logKernelOut
  );
}


// Cumulative-probit likelihood with a separate, equally-sized quadrature rule
// for every response pattern. Node-specific means and quadrature vectors are
// packed in consecutive pattern blocks.
// [[Rcpp::export]]
double adaptiveOrdinalLogLikCpp(const Rcpp::IntegerMatrix& codes,
                                 const arma::mat& mu,
                                 const Rcpp::List& thresholds,
                                 const arma::vec& quadWeights,
                                 const arma::vec& logAdjust,
                                 const arma::vec& patternWeights,
                                 const int nodesPerPattern) {
  const arma::uword n = codes.nrow(), jn = codes.ncol();
  const arma::uword qn = static_cast<arma::uword>(nodesPerPattern);
  if (nodesPerPattern < 1 || mu.n_rows != n * qn || mu.n_cols != jn)
    Rcpp::stop("Packed adaptive ordinal means have incompatible dimensions.");
  if (quadWeights.n_elem != n * qn || logAdjust.n_elem != n * qn ||
      patternWeights.n_elem != n || thresholds.size() != static_cast<R_xlen_t>(jn))
    Rcpp::stop("Packed adaptive ordinal inputs have incompatible dimensions.");
  std::vector<arma::vec> tau(jn);
  for (arma::uword j = 0; j < jn; ++j) tau[j] = Rcpp::as<arma::vec>(thresholds[j]);
  double logLik = 0.0;
  arma::vec logKernel(qn), mass(qn);
  for (arma::uword i = 0; i < n; ++i) {
    logKernel.zeros();
    const arma::uword base = i * qn;
    for (arma::uword j = 0; j < jn; ++j) {
      const int code = codes(i, j);
      if (code == NA_INTEGER || code <= 0) continue;
      const int categories = static_cast<int>(tau[j].n_elem) + 1;
      if (code > categories) Rcpp::stop("Category code exceeds its threshold count.");
      const double lo = code == 1 ? R_NegInf : tau[j][code - 2];
      const double hi = code == categories ? R_PosInf : tau[j][code - 1];
      for (arma::uword q = 0; q < qn; ++q) {
        const double eta = mu(base + q, j);
        const double Flo = std::isinf(lo) ? 0.0 : normal_cdf(lo - eta);
        const double Fhi = std::isinf(hi) ? 1.0 : normal_cdf(hi - eta);
        logKernel[q] += std::log(std::max(Fhi - Flo, DBL_MIN));
      }
    }
    double scale = R_NegInf;
    for (arma::uword q = 0; q < qn; ++q)
      if (quadWeights[base + q] != 0.0)
        scale = std::max(scale, logKernel[q] + logAdjust[base + q] +
                         std::log(std::abs(quadWeights[base + q])));
    double integral = 0.0;
    for (arma::uword q = 0; q < qn; ++q)
      integral += quadWeights[base + q] *
        std::exp(logKernel[q] + logAdjust[base + q] - scale);
    if (!(integral > 0.0) || !std::isfinite(integral)) return R_NegInf;
    logLik += patternWeights[i] * (scale + std::log(integral));
  }
  return logLik;
}


// Fixed-quadrature kernel for mixed cumulative-probit and Gaussian indicators.
// Ordinal and continuous indicators may both contain missing values.  The
// continuous covariance is node-specific; its score uses the same ordering as
// `covariance` and is returned as one matrix per quadrature node.
// [[Rcpp::export]]
Rcpp::List mixedLogLikCpp(const Rcpp::IntegerMatrix& codes,
                          const arma::mat& categoricalMu,
                          const Rcpp::List& thresholds,
                          const arma::mat& continuousData,
                          const arma::mat& continuousMu,
                          const Rcpp::List& covariance,
                          const arma::vec& quadWeights,
                          const arma::vec& logAdjust,
                          const arma::vec& rowWeights,
                          const bool gradient = false,
                          const bool details = false,
                          const bool posteriorDetails = true) {
  const arma::uword n = continuousData.n_rows;
  const arma::uword jc = codes.ncol();
  const arma::uword jx = continuousData.n_cols;
  const arma::uword qn = categoricalMu.n_rows;
  if (static_cast<arma::uword>(codes.nrow()) != n)
    Rcpp::stop("Categorical and continuous data must have equal rows.");
  if (categoricalMu.n_cols != jc || continuousMu.n_rows != qn ||
      continuousMu.n_cols != jx)
    Rcpp::stop("Node-specific means have incompatible dimensions.");
  if (thresholds.size() != static_cast<R_xlen_t>(jc))
    Rcpp::stop("One threshold vector is required per categorical indicator.");
  if (covariance.size() != static_cast<R_xlen_t>(qn))
    Rcpp::stop("One continuous covariance matrix is required per node.");
  if (quadWeights.n_elem != qn || logAdjust.n_elem != qn)
    Rcpp::stop("Quadrature vectors must match the node count.");
  if (rowWeights.n_elem != n)
    Rcpp::stop("Row weights must match the data rows.");

  std::vector<arma::vec> tau(jc);
  arma::uvec offsets(jc + 1, arma::fill::zeros);
  for (arma::uword j = 0; j < jc; ++j) {
    tau[j] = Rcpp::as<arma::vec>(thresholds[j]);
    offsets[j + 1] = offsets[j] + tau[j].n_elem;
  }
  std::vector<arma::mat> sigma(qn);
  for (arma::uword q = 0; q < qn; ++q) {
    sigma[q] = Rcpp::as<arma::mat>(covariance[q]);
    if (sigma[q].n_rows != jx || sigma[q].n_cols != jx)
      Rcpp::stop("Continuous covariance has incompatible dimensions.");
  }

  arma::mat categoricalScore(qn, jc, arma::fill::zeros);
  arma::mat continuousScore(qn, jx, arma::fill::zeros);
  std::vector<arma::mat> covarianceScore(qn, arma::mat(jx, jx, arma::fill::zeros));
  arma::vec thresholdScore(offsets[jc], arma::fill::zeros);
  arma::mat dcat(qn, jc, arma::fill::zeros);
  arma::mat upperDerivative(qn, jc, arma::fill::zeros);
  arma::mat lowerDerivative(qn, jc, arma::fill::zeros);
  std::vector<arma::vec> dcontinuous(qn, arma::vec(jx, arma::fill::zeros));
  std::vector<arma::mat> dcovariance(qn, arma::mat(jx, jx, arma::fill::zeros));
  arma::ivec category(jc, arma::fill::zeros);
  arma::vec logKernel(qn, arma::fill::zeros);
  arma::vec signedMass(qn, arma::fill::zeros);
  arma::mat posteriorOut(n, qn, arma::fill::zeros);
  arma::mat logKernelOut(n, qn, arma::fill::zeros);
  double logLik = 0.0;

  // Covariance submatrices depend only on the node and the continuous
  // missingness pattern, not on the observed values. Factor each unique
  // pattern once instead of once per row and node.
  std::unordered_map<std::string, arma::uword> patternMap;
  std::vector<arma::uvec> patternIndex;
  arma::uvec rowPattern(n);
  for (arma::uword i = 0; i < n; ++i) {
    std::string key;
    key.reserve(jx);
    std::vector<arma::uword> observed;
    for (arma::uword j = 0; j < jx; ++j) {
      const bool present = std::isfinite(continuousData(i, j));
      key.push_back(present ? '1' : '0');
      if (present) observed.push_back(j);
    }
    const auto found = patternMap.find(key);
    if (found != patternMap.end()) rowPattern[i] = found->second;
    else {
      const arma::uword id = patternIndex.size();
      patternMap[key] = id;
      rowPattern[i] = id;
      arma::uvec idx(observed.size());
      for (arma::uword h = 0; h < observed.size(); ++h) idx[h] = observed[h];
      patternIndex.push_back(idx);
    }
  }
  std::vector<std::vector<arma::mat> > patternL(
    patternIndex.size(), std::vector<arma::mat>(qn));
  std::vector<std::vector<arma::mat> > patternInv(
    patternIndex.size(), std::vector<arma::mat>(qn));
  std::vector<std::vector<double> > patternLogDet(
    patternIndex.size(), std::vector<double>(qn, 0.0));
  for (arma::uword h = 0; h < patternIndex.size(); ++h) {
    const arma::uvec& idx = patternIndex[h];
    if (idx.n_elem == 0) continue;
    for (arma::uword q = 0; q < qn; ++q) {
      const arma::mat S = sigma[q].submat(idx, idx);
      if (!arma::chol(patternL[h][q], S, "lower"))
        return Rcpp::List::create(Rcpp::Named("logLik") = R_NegInf);
      patternLogDet[h][q] = 2.0 * arma::accu(arma::log(patternL[h][q].diag()));
      if (gradient)
        patternInv[h][q] = arma::solve(arma::trimatu(patternL[h][q].t()),
          arma::solve(arma::trimatl(patternL[h][q]),
                      arma::eye(idx.n_elem, idx.n_elem)));
    }
  }

  for (arma::uword i = 0; i < n; ++i) {
    logKernel.zeros();
    if (gradient) {
      dcat.zeros(); upperDerivative.zeros(); lowerDerivative.zeros();
      for (arma::uword q = 0; q < qn; ++q) {
        dcontinuous[q].zeros(); dcovariance[q].zeros();
      }
    }
    for (arma::uword j = 0; j < jc; ++j) {
      const int code = codes(i, j);
      category[j] = code;
      if (code == NA_INTEGER || code <= 0) continue;
      const int categories = static_cast<int>(tau[j].n_elem) + 1;
      if (code > categories) Rcpp::stop("Category code exceeds its threshold count.");
      const double lo = code == 1 ? R_NegInf : tau[j][code - 2];
      const double hi = code == categories ? R_PosInf : tau[j][code - 1];
      for (arma::uword q = 0; q < qn; ++q) {
        const double eta = categoricalMu(q, j);
        const double Flo = std::isinf(lo) ? 0.0 : normal_cdf(lo - eta);
        const double Fhi = std::isinf(hi) ? 1.0 : normal_cdf(hi - eta);
        const double probability = std::max(Fhi - Flo, DBL_MIN);
        logKernel[q] += std::log(probability);
        if (gradient) {
          const double flo = normal_density(lo - eta);
          const double fhi = normal_density(hi - eta);
          dcat(q, j) = (flo - fhi) / probability;
          lowerDerivative(q, j) = -flo / probability;
          upperDerivative(q, j) = fhi / probability;
        }
      }
    }

    const arma::uword pattern = rowPattern[i];
    const arma::uvec& idx = patternIndex[pattern];
    if (idx.n_elem > 0) {
      arma::vec y(idx.n_elem);
      for (arma::uword h = 0; h < idx.n_elem; ++h)
        y[h] = continuousData(i, idx[h]);
      for (arma::uword q = 0; q < qn; ++q) {
        const arma::mat& L = patternL[pattern][q];
        arma::vec muObserved(idx.n_elem);
        for (arma::uword h = 0; h < idx.n_elem; ++h)
          muObserved[h] = continuousMu(q, idx[h]);
        const arma::vec diff = y - muObserved;
        const arma::vec z = arma::solve(arma::trimatl(L), diff);
        logKernel[q] -= 0.5 * (idx.n_elem * std::log(2.0 * M_PI) +
          patternLogDet[pattern][q] + arma::dot(z, z));
        if (gradient) {
          const arma::vec invDiff = arma::solve(arma::trimatu(L.t()), z);
          dcontinuous[q].elem(idx) = invDiff;
          dcovariance[q].submat(idx, idx) =
            0.5 * (invDiff * invDiff.t() - patternInv[pattern][q]);
        }
      }
    }

    double scale = R_NegInf;
    for (arma::uword q = 0; q < qn; ++q) if (quadWeights[q] != 0.0)
      scale = std::max(scale, logKernel[q] + logAdjust[q] +
                       std::log(std::abs(quadWeights[q])));
    if (!std::isfinite(scale))
      return Rcpp::List::create(Rcpp::Named("logLik") = R_NegInf);
    double integralScaled = 0.0;
    for (arma::uword q = 0; q < qn; ++q) {
      signedMass[q] = quadWeights[q] * std::exp(logKernel[q] + logAdjust[q] - scale);
      integralScaled += signedMass[q];
    }
    if (!(integralScaled > 0.0) || !std::isfinite(integralScaled))
      return Rcpp::List::create(Rcpp::Named("logLik") = R_NegInf);
    const double frequency = rowWeights[i];
    logLik += frequency * (scale + std::log(integralScaled));
    if (gradient || details) {
      logKernelOut.row(i) = logKernel.t();
      posteriorOut.row(i) = (signedMass / integralScaled).t();
    }
    if (!gradient) continue;

    for (arma::uword q = 0; q < qn; ++q) {
      const double posterior = frequency * signedMass[q] / integralScaled;
      categoricalScore.row(q) += posterior * dcat.row(q);
      continuousScore.row(q) += posterior * dcontinuous[q].t();
      covarianceScore[q] += posterior * dcovariance[q];
      for (arma::uword j = 0; j < jc; ++j) {
        const int code = category[j];
        if (code == NA_INTEGER || code <= 0) continue;
        const int categories = static_cast<int>(tau[j].n_elem) + 1;
        if (code > 1)
          thresholdScore[offsets[j] + code - 2] += posterior * lowerDerivative(q, j);
        if (code < categories)
          thresholdScore[offsets[j] + code - 1] += posterior * upperDerivative(q, j);
      }
    }
  }
  Rcpp::List covarianceScoreOut(qn);
  for (arma::uword q = 0; q < qn; ++q) covarianceScoreOut[q] = covarianceScore[q];
  if (!gradient && !details)
    return Rcpp::List::create(Rcpp::Named("logLik") = logLik);
  if (!gradient)
    return posteriorDetails ? Rcpp::List::create(
        Rcpp::Named("logLik") = logLik,
        Rcpp::Named("posterior") = posteriorOut,
        Rcpp::Named("log.kernel") = logKernelOut
      ) : Rcpp::List::create(
        Rcpp::Named("logLik") = logLik,
        Rcpp::Named("log.kernel") = logKernelOut
      );
  return Rcpp::List::create(
    Rcpp::Named("logLik") = logLik,
    Rcpp::Named("categorical.score") = categoricalScore,
    Rcpp::Named("continuous.score") = continuousScore,
    Rcpp::Named("covariance.score") = covarianceScoreOut,
    Rcpp::Named("threshold.score") = thresholdScore,
    Rcpp::Named("posterior") = posteriorOut,
    Rcpp::Named("log.kernel") = logKernelOut
  );
}


// Mixed ordinal/Gaussian likelihood with one quadrature rule per row. All
// node-dependent quantities use consecutive row blocks of `nodesPerRow`.
// [[Rcpp::export]]
double adaptiveMixedLogLikCpp(const Rcpp::IntegerMatrix& codes,
                               const arma::mat& categoricalMu,
                               const Rcpp::List& thresholds,
                               const arma::mat& continuousData,
                               const arma::mat& continuousMu,
                               const Rcpp::List& covariance,
                               const arma::vec& quadWeights,
                               const arma::vec& logAdjust,
                               const arma::vec& rowWeights,
                               const int nodesPerRow) {
  const arma::uword n = continuousData.n_rows;
  const arma::uword jc = codes.ncol(), jx = continuousData.n_cols;
  const arma::uword qn = static_cast<arma::uword>(nodesPerRow);
  const arma::uword packed = n * qn;
  if (nodesPerRow < 1 || static_cast<arma::uword>(codes.nrow()) != n ||
      categoricalMu.n_rows != packed || categoricalMu.n_cols != jc ||
      continuousMu.n_rows != packed || continuousMu.n_cols != jx ||
      covariance.size() != static_cast<R_xlen_t>(packed) ||
      quadWeights.n_elem != packed || logAdjust.n_elem != packed ||
      rowWeights.n_elem != n || thresholds.size() != static_cast<R_xlen_t>(jc))
    Rcpp::stop("Packed adaptive mixed-response inputs have incompatible dimensions.");
  std::vector<arma::vec> tau(jc);
  for (arma::uword j = 0; j < jc; ++j) tau[j] = Rcpp::as<arma::vec>(thresholds[j]);
  double logLik = 0.0;
  arma::vec logKernel(qn);
  for (arma::uword i = 0; i < n; ++i) {
    logKernel.zeros();
    const arma::uword base = i * qn;
    for (arma::uword j = 0; j < jc; ++j) {
      const int code = codes(i, j);
      if (code == NA_INTEGER || code <= 0) continue;
      const int categories = static_cast<int>(tau[j].n_elem) + 1;
      if (code > categories) Rcpp::stop("Category code exceeds its threshold count.");
      const double lo = code == 1 ? R_NegInf : tau[j][code - 2];
      const double hi = code == categories ? R_PosInf : tau[j][code - 1];
      for (arma::uword q = 0; q < qn; ++q) {
        const double eta = categoricalMu(base + q, j);
        const double Flo = std::isinf(lo) ? 0.0 : normal_cdf(lo - eta);
        const double Fhi = std::isinf(hi) ? 1.0 : normal_cdf(hi - eta);
        logKernel[q] += std::log(std::max(Fhi - Flo, DBL_MIN));
      }
    }
    std::vector<arma::uword> observed;
    for (arma::uword j = 0; j < jx; ++j)
      if (std::isfinite(continuousData(i, j))) observed.push_back(j);
    arma::uvec idx(observed.size());
    arma::vec y(observed.size());
    for (arma::uword h = 0; h < observed.size(); ++h) {
      idx[h] = observed[h]; y[h] = continuousData(i, observed[h]);
    }
    if (idx.n_elem > 0) for (arma::uword q = 0; q < qn; ++q) {
      const arma::mat sigma = Rcpp::as<arma::mat>(covariance[base + q]);
      const arma::mat S = sigma.submat(idx, idx);
      arma::mat L;
      if (!arma::chol(L, S, "lower")) return R_NegInf;
      arma::vec mean(idx.n_elem);
      for (arma::uword h = 0; h < idx.n_elem; ++h)
        mean[h] = continuousMu(base + q, idx[h]);
      const arma::vec z = arma::solve(arma::trimatl(L), y - mean);
      logKernel[q] -= 0.5 * (idx.n_elem * std::log(2.0 * M_PI) +
        2.0 * arma::accu(arma::log(L.diag())) + arma::dot(z, z));
    }
    double scale = R_NegInf;
    for (arma::uword q = 0; q < qn; ++q)
      if (quadWeights[base + q] != 0.0)
        scale = std::max(scale, logKernel[q] + logAdjust[base + q] +
                         std::log(std::abs(quadWeights[base + q])));
    double integral = 0.0;
    for (arma::uword q = 0; q < qn; ++q)
      integral += quadWeights[base + q] *
        std::exp(logKernel[q] + logAdjust[base + q] - scale);
    if (!(integral > 0.0) || !std::isfinite(integral)) return R_NegInf;
    logLik += rowWeights[i] * (scale + std::log(integral));
  }
  return logLik;
}


// Node-specific conditional latent means used by the categorical kernel.
// `explicit` contains one-based positions in the full independent innovation
// vector. Analytically integrated innovations are fixed at their zero mean.
// [[Rcpp::export]]
arma::mat lmsCatLatentMeansCpp(const Rcpp::List& matrices,
                               const arma::mat& nodes,
                               const Rcpp::IntegerVector& explicitIndex,
                               const int fullDimension) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat covZetaXi = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat omegaXiXi = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat omegaEtaXi = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);

  const arma::uword p = A.n_rows;
  const arma::uword r = psi.n_rows;
  if (nodes.n_cols != static_cast<arma::uword>(explicitIndex.size()))
    Rcpp::stop("Nodes do not match the explicit integration dimensions.");
  if (fullDimension != static_cast<int>(p + r))
    Rcpp::stop("Full integration dimension does not match the model.");

  arma::mat zproj(r, p, arma::fill::zeros);
  arma::mat innovationRoot(r, r, arma::fill::zeros);
  if (r > 0) {
    zproj = arma::solve(A, covZetaXi.t()).t();
    arma::mat orthogonalPsi = psi - zproj * zproj.t();
    orthogonalPsi = 0.5 * (orthogonalPsi + orthogonalPsi.t());
    if (!arma::chol(innovationRoot, orthogonalPsi, "lower"))
      Rcpp::stop("The conditional latent covariance is not positive definite.");
  }

  arma::mat result(nodes.n_rows, p + r, arma::fill::zeros);
  const arma::mat identityEta = arma::eye(r, r);
  for (arma::uword q = 0; q < nodes.n_rows; ++q) {
    arma::vec u(fullDimension, arma::fill::zeros);
    for (R_xlen_t k = 0; k < explicitIndex.size(); ++k) {
      const int position = explicitIndex[k] - 1;
      if (position < 0 || position >= fullDimension)
        Rcpp::stop("Invalid explicit integration index.");
      u[position] = nodes(q, k);
    }
    const arma::vec ux = u.head(p);
    const arma::vec xi = beta0 + A * ux;
    result(q, arma::span(0, p - 1)) = xi.t();
    if (r == 0) continue;

    const arma::vec ue = u.tail(r);
    const arma::vec zeta = zproj * ux + innovationRoot * ue;
    const arma::mat x(xi);
    const arma::mat kz = arma::kron(identityEta, x);
    const arma::mat B = identityEta - gammaEta - kz.t() * omegaEtaXi;
    const arma::vec rhs = alpha + gammaXi * xi +
      kz.t() * omegaXiXi * xi + zeta;
    const arma::vec eta = arma::solve(B, rhs);
    result(q, arma::span(p, p + r - 1)) = eta.t();
  }
  return result;
}


// Conditional latent means and the covariance contributed by analytically
// integrated (linear) innovations. This is the compiled equivalent of
// `lmsCatReducedMoments` and accepts an arbitrary packed node matrix.
// [[Rcpp::export]]
Rcpp::List lmsCatReducedMomentsCpp(const Rcpp::List& matrices,
                                   const arma::mat& nodes,
                                   const Rcpp::IntegerVector& explicitIndex,
                                   const Rcpp::IntegerVector& analyticIndex,
                                   const int fullDimension) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat covZetaXi = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat gammaXi = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat gammaEta = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat omegaXiXi = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat omegaEtaXi = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);
  const arma::uword p = A.n_rows, r = psi.n_rows, latentDimension = p + r;
  if (nodes.n_cols != static_cast<arma::uword>(explicitIndex.size()) ||
      fullDimension != static_cast<int>(latentDimension))
    Rcpp::stop("Reduced-moment integration dimensions are incompatible.");
  arma::mat zproj(r, p, arma::fill::zeros), innovationRoot(r, r, arma::fill::zeros);
  if (r > 0) {
    zproj = arma::solve(A, covZetaXi.t()).t();
    arma::mat orthogonalPsi = psi - zproj * zproj.t();
    orthogonalPsi = 0.5 * (orthogonalPsi + orthogonalPsi.t());
    if (!arma::chol(innovationRoot, orthogonalPsi, "lower"))
      Rcpp::stop("The conditional latent covariance is not positive definite.");
  }
  const arma::mat identityEta = arma::eye(r, r);
  const auto stateOne = [&](const arma::vec& u) -> arma::vec {
    const arma::vec ux = u.head(p);
    const arma::vec xi = beta0 + A * ux;
    if (r == 0) return xi;
    const arma::vec ue = u.tail(r);
    const arma::vec zeta = zproj * ux + innovationRoot * ue;
    const arma::mat x(xi);
    const arma::mat kz = arma::kron(identityEta, x);
    const arma::mat B = identityEta - gammaEta - kz.t() * omegaEtaXi;
    const arma::vec rhs = alpha + gammaXi * xi +
      kz.t() * omegaXiXi * xi + zeta;
    return arma::join_cols(xi, arma::solve(B, rhs));
  };
  arma::mat means(nodes.n_rows, latentDimension, arma::fill::zeros);
  Rcpp::List covariance(nodes.n_rows);
  for (arma::uword q = 0; q < nodes.n_rows; ++q) {
    arma::vec u(fullDimension, arma::fill::zeros);
    for (R_xlen_t h = 0; h < explicitIndex.size(); ++h)
      u[explicitIndex[h] - 1] = nodes(q, h);
    const arma::vec base = stateOne(u);
    means.row(q) = base.t();
    arma::mat J(latentDimension, analyticIndex.size(), arma::fill::zeros);
    for (R_xlen_t h = 0; h < analyticIndex.size(); ++h) {
      arma::vec shifted = u;
      shifted[analyticIndex[h] - 1] = 1.0;
      J.col(h) = stateOne(shifted) - base;
    }
    covariance[q] = J * J.t();
  }
  return Rcpp::List::create(Rcpp::Named("mean") = means,
                            Rcpp::Named("cov") = covariance);
}


// Analytical forward sensitivities of the categorical linear predictors with
// respect to model-matrix entries. The returned score is on the matrix-entry
// scale; R applies the model constraint Jacobian to obtain free-parameter
// derivatives.
// [[Rcpp::export]]
arma::vec lmsCatMatrixScoreCpp(const Rcpp::List& matrices,
                               const arma::mat& nodes,
                               const Rcpp::IntegerVector& explicitIndex,
                               const int fullDimension,
                               const arma::mat& predictorScore,
                               const Rcpp::List& covarianceScore,
                               const Rcpp::IntegerVector& continuousIndex,
                               const Rcpp::IntegerVector& analyticIndex,
                               const arma::uvec& block,
                               const arma::uvec& row,
                               const arma::uvec& col,
                               const arma::uvec& symmetric) {
  const arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  const arma::mat psi = Rcpp::as<arma::mat>(matrices["psi"]);
  const arma::mat cov = Rcpp::as<arma::mat>(matrices["covZetaXi"]);
  const arma::mat lambda = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  const arma::mat theta = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  const arma::vec tau = Rcpp::as<arma::vec>(matrices["tauX"]);
  const arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  const arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  const arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  const arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  const arma::vec beta0 = Rcpp::as<arma::vec>(matrices["beta0"]);
  const arma::vec alpha = Rcpp::as<arma::vec>(matrices["alpha"]);
  const arma::uword p = A.n_rows;
  const arma::uword r = psi.n_rows;
  const arma::uword jobs = lambda.n_rows;
  const arma::uword npar = block.n_elem;
  if (predictorScore.n_rows != nodes.n_rows || predictorScore.n_cols != jobs)
    Rcpp::stop("Predictor score has incompatible dimensions.");
  if (row.n_elem != npar || col.n_elem != npar || symmetric.n_elem != npar)
    Rcpp::stop("Parameter-location vectors have incompatible lengths.");
  const bool hasCovarianceScore = covarianceScore.size() > 0;
  if (hasCovarianceScore && covarianceScore.size() !=
      static_cast<R_xlen_t>(nodes.n_rows))
    Rcpp::stop("Covariance score must contain one matrix per node.");
  arma::uvec continuous(continuousIndex.size());
  arma::uvec analytic(analyticIndex.size());
  for (R_xlen_t h = 0; h < continuousIndex.size(); ++h) {
    if (continuousIndex[h] <= 0 || continuousIndex[h] > static_cast<int>(jobs))
      Rcpp::stop("Invalid continuous-indicator index.");
    continuous[h] = continuousIndex[h] - 1;
  }
  for (R_xlen_t h = 0; h < analyticIndex.size(); ++h) {
    if (analyticIndex[h] <= 0 || analyticIndex[h] > fullDimension)
      Rcpp::stop("Invalid analytic integration index.");
    analytic[h] = analyticIndex[h] - 1;
  }

  const arma::mat invAt = arma::inv(A).t();
  const arma::mat Z = cov * invAt;
  arma::mat S = psi - Z * Z.t();
  S = 0.5 * (S + S.t());
  arma::mat L;
  if (r > 0 && !arma::chol(L, S, "lower"))
    Rcpp::stop("The conditional latent covariance is not positive definite.");
  const arma::mat Ie = arma::eye(r, r);
  arma::vec result(npar, arma::fill::zeros);

  for (arma::uword k = 0; k < npar; ++k) {
    arma::mat dA(A.n_rows, A.n_cols, arma::fill::zeros);
    arma::mat dPsi(psi.n_rows, psi.n_cols, arma::fill::zeros);
    arma::mat dCov(cov.n_rows, cov.n_cols, arma::fill::zeros);
    arma::mat dLambda(lambda.n_rows, lambda.n_cols, arma::fill::zeros);
    arma::mat dTheta(theta.n_rows, theta.n_cols, arma::fill::zeros);
    arma::vec dTau(tau.n_elem, arma::fill::zeros);
    arma::mat dGx(Gx.n_rows, Gx.n_cols, arma::fill::zeros);
    arma::mat dGe(Ge.n_rows, Ge.n_cols, arma::fill::zeros);
    arma::mat dOxx(Oxx.n_rows, Oxx.n_cols, arma::fill::zeros);
    arma::mat dOex(Oex.n_rows, Oex.n_cols, arma::fill::zeros);
    arma::vec dBeta(beta0.n_elem, arma::fill::zeros);
    arma::vec dAlpha(alpha.n_elem, arma::fill::zeros);

    const arma::uword rr = row[k], cc = col[k];
    switch (block[k]) {
    case 0: dLambda(rr, cc) = 1.0; break;
    case 2: dTau[rr] = 1.0; break;
    case 4:
      dTheta(rr, cc) = 1.0;
      if (symmetric[k] && rr != cc) dTheta(cc, rr) = 1.0;
      break;
    case 6: dA(rr, cc) = 1.0; break;
    case 7:
      dPsi(rr, cc) = 1.0;
      if (symmetric[k] && rr != cc) dPsi(cc, rr) = 1.0;
      break;
    case 8: dAlpha[rr] = 1.0; break;
    case 9: dBeta[rr] = 1.0; break;
    case 10: dGx(rr, cc) = 1.0; break;
    case 11: dGe(rr, cc) = 1.0; break;
    case 12:
      dOxx(rr, cc) = 1.0;
      if (symmetric[k] && rr != cc) dOxx(cc, rr) = 1.0;
      break;
    case 13: dOex(rr, cc) = 1.0; break;
    case 17: dCov(rr, cc) = 1.0; break;
    default: continue; // residual and composite-only blocks do not enter here
    }

    const arma::mat dZ = dCov * invAt - Z * dA.t() * invAt;
    arma::mat dL(r, r, arma::fill::zeros);
    if (r > 0 && (arma::accu(arma::abs(dPsi)) > 0.0 ||
                  arma::accu(arma::abs(dZ)) > 0.0)) {
      arma::mat dS = dPsi - dZ * Z.t() - Z * dZ.t();
      arma::mat E = arma::solve(arma::trimatl(L), dS);
      E = arma::solve(arma::trimatl(L), E.t()).t();
      arma::mat Phi = arma::trimatl(E);
      Phi.diag() *= 0.5;
      dL = L * Phi;
    }

    const auto stateSensitivity = [&](const arma::vec& u) {
      const arma::vec ux = u.head(p);
      const arma::vec xi = beta0 + A * ux;
      const arma::vec dxi = dBeta + dA * ux;
      arma::vec latent = xi;
      arma::vec dlatent = dxi;
      if (r > 0) {
        const arma::vec ue = u.tail(r);
        const arma::vec zeta = Z * ux + L * ue;
        const arma::vec dzeta = dZ * ux + dL * ue;
        const arma::mat x(xi), dx(dxi);
        const arma::mat K = arma::kron(Ie, x);
        const arma::mat dK = arma::kron(Ie, dx);
        const arma::mat B = Ie - Ge - K.t() * Oex;
        const arma::mat dB = -dGe - dK.t() * Oex - K.t() * dOex;
        const arma::vec rhs = alpha + Gx * xi + K.t() * Oxx * xi + zeta;
        const arma::vec eta = arma::solve(B, rhs);
        const arma::vec drhs = dAlpha + dGx * xi + Gx * dxi +
          dK.t() * Oxx * xi + K.t() * dOxx * xi +
          K.t() * Oxx * dxi + dzeta;
        const arma::vec deta = arma::solve(B, drhs - dB * eta);
        latent = arma::join_cols(xi, eta);
        dlatent = arma::join_cols(dxi, deta);
      }
      return std::make_pair(latent, dlatent);
    };

    for (arma::uword q = 0; q < nodes.n_rows; ++q) {
      arma::vec u(fullDimension, arma::fill::zeros);
      for (R_xlen_t h = 0; h < explicitIndex.size(); ++h)
        u[explicitIndex[h] - 1] = nodes(q, h);
      const std::pair<arma::vec, arma::vec> base = stateSensitivity(u);
      const arma::vec& latent = base.first;
      const arma::vec& dlatent = base.second;
      const arma::vec dmu = dTau + dLambda * latent + lambda * dlatent;
      result[k] += arma::dot(predictorScore.row(q), dmu.t());

      if (hasCovarianceScore && continuous.n_elem > 0) {
        arma::mat J(latent.n_elem, analytic.n_elem, arma::fill::zeros);
        arma::mat dJ(latent.n_elem, analytic.n_elem, arma::fill::zeros);
        for (arma::uword h = 0; h < analytic.n_elem; ++h) {
          arma::vec shifted = u;
          shifted[analytic[h]] = 1.0;
          const std::pair<arma::vec, arma::vec> shiftedState =
            stateSensitivity(shifted);
          J.col(h) = shiftedState.first - latent;
          dJ.col(h) = shiftedState.second - dlatent;
        }
        const arma::mat latentCov = J * J.t();
        const arma::mat dLatentCov = dJ * J.t() + J * dJ.t();
        const arma::mat dObservedCov =
          dLambda * latentCov * lambda.t() +
          lambda * dLatentCov * lambda.t() +
          lambda * latentCov * dLambda.t() + dTheta;
        const arma::mat dContinuousCov = dObservedCov.submat(continuous, continuous);
        const arma::mat score = Rcpp::as<arma::mat>(covarianceScore[q]);
        if (score.n_rows != continuous.n_elem || score.n_cols != continuous.n_elem)
          Rcpp::stop("Covariance score has incompatible dimensions.");
        result[k] += arma::accu(score % dContinuousCov);
      }
    }
  }
  return result;
}
