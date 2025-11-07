#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>
#include <math.h>

#include "utils.h"
#include "mvnorm.h"
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

arma::vec dnorm(const arma::vec x, const double mu, const double sd, const bool log = false) {
  const arma::vec e = x - mu;
  const arma::vec ldens = (- 0.5) * std::log(2 * M_PI) - std::log(sd) - (e % e)/(2 * sd * sd);
  return log ? ldens: arma::exp(ldens);
}

inline int resolve_ncores(const int requested) {
#ifdef _OPENMP
  if (requested <= 0)
    return std::max(1, omp_get_max_threads());
  return std::max(1, std::min(requested, omp_get_max_threads()));
#else
  (void)requested;
  return 1;
#endif
}


arma::vec pnormOrderedProbit(const arma::vec v, // Expected value for latent response variable
                             const arma::vec lower, // lower thresholds for each response
                             const arma::vec upper, // upper thresholds for each response
                             const double mean = 0,
                             const double sd = 1,
                             const bool log = false)
{
  const arma::uword n = v.n_elem;

  if (lower.n_elem != n || upper.n_elem != n)
    Rcpp::stop("`lower` and `upper` must be the same length as `v`.");

  if (sd <= 0)
    Rcpp::stop("Standard deviation `sd` must be positive.");

  arma::vec out(n);
  const double sqrt2 = std::sqrt(2.0);
  const double invSdSqrt2 = 1.0 / (sd * sqrt2);
  const arma::vec mu = mean + v;

  arma::vec zLower = (lower - mu) * invSdSqrt2;
  arma::vec zUpper = (upper - mu) * invSdSqrt2;

  arma::vec lowerCdf = 0.5 * arma::erfc(-zLower);
  arma::vec upperCdf = 0.5 * arma::erfc(-zUpper);

  arma::vec prob = upperCdf - lowerCdf;
  arma::uvec nonPositive = arma::find(prob <= 0.0);
  prob.elem(nonPositive).fill(arma::datum::nan);

  return log ? arma::log(prob) : prob;
}


// Gsem model for a single group/cluster
struct GSEM_ModelGroup {
  arma::mat Lambda, tau, Thresholds, Ie, Gamma, Psi, alpha, Theta, W;

  std::vector<arma::mat> Z, Y;
  std::vector<arma::uvec> colIdxPatterns;
  
  arma::uvec n;

  struct FreeEntry {
    arma::uword row = 0;
    arma::uword col = 0;
    bool symmetric = false;
    std::string label;
  };

  arma::mat cholPsi, B, Bt, LambdaT;
  arma::rowvec alphaRow;
  arma::vec thetaStd, thetaInvVar, thetaLogNorm;
  std::vector<arma::vec> thresholdRows;
  bool cholPsiValid = true;
  std::vector<FreeEntry> freeLambda, freeTau, freeGamma, freeAlpha,
    freePsi, freeTheta, freeThresholds;
  std::vector<std::string> freeLabelOrder;

  struct PredictorCacheEntry {
    arma::mat latent;
    arma::mat S;
    arma::mat V;
    bool valid = false;
  };

  mutable std::vector<PredictorCacheEntry> predictorCache;
  mutable bool predictorCacheInitialized = false;

  unsigned k = 0, p = 0, N = 0;
  arma::uvec isordered; // if 0 variable is not ordered, if > 0, then it's the 1 based index
                        // denoting the row in `thresholds`

  explicit GSEM_ModelGroup(const Rcpp::List& modFilled) {

    Rcpp::List matrices  = modFilled["matrices"];
    Rcpp::List info      = modFilled["info"];

    Rcpp::List quad      = modFilled["quad"];
    Rcpp::List quad_i    = quad["n"];

    Rcpp::List dataR     = modFilled["data"];
    Rcpp::List dataSplit = dataR["data.split"];
    Rcpp::List colidxR   = dataR["colidx"];

    n = Rcpp::as<arma::uvec>(dataR["n"]); // Observations per patter
    N = arma::sum(n);
    k = Rcpp::as<unsigned>(dataR["k"]); // Cols
    p = Rcpp::as<unsigned>(dataR["p"]); // patterns

    // one-liners, no loops
    Ie         = Rcpp::as<arma::mat>(matrices["Ieta"]);
    Lambda     = Rcpp::as<arma::mat>(matrices["lambda"]);
    tau        = Rcpp::as<arma::mat>(matrices["tau"]);
    Gamma      = Rcpp::as<arma::mat>(matrices["gamma"]);
    alpha      = Rcpp::as<arma::mat>(matrices["alpha"]);
    Psi        = Rcpp::as<arma::mat>(matrices["psi"]);
    Theta      = Rcpp::as<arma::mat>(matrices["theta"]);
    Thresholds = Rcpp::as<arma::mat>(matrices["thresholds"]);
    isordered  = Rcpp::as<arma::uvec>(matrices["isordered"]);

    W = Rcpp::as<arma::mat>(quad["w"]);
    Z = std::vector<arma::mat>(quad_i.length());

    for (int q = 0; q < quad_i.length(); q++)
      Z[q] = Rcpp::as<arma::mat>(quad_i[q]);

    Y              = std::vector<arma::mat>(p);
    colIdxPatterns = std::vector<arma::uvec>(p);

    for (int pi = 0; pi < p; pi++) {
      Y[pi]              = Rcpp::as<arma::mat>(dataSplit[pi]);
      colIdxPatterns[pi] = Rcpp::as<arma::uvec>(colidxR[pi]) - 1L; // Make zero based
    }

    initializeCaches();
    initializeFreeParams(modFilled);
  }
  
  explicit GSEM_ModelGroup() { // Empty initilizer
    k = 0, N = 0, p = 0;
  }

  void initializeCaches() {
    const arma::mat IeMinusGamma = Ie - Gamma;
    const arma::uword latentDim = IeMinusGamma.n_rows;

    cholPsiValid = true;
    if (Psi.n_rows > 0) {
      cholPsiValid = arma::chol(cholPsi, Psi, "upper");
      if (!cholPsiValid)
        cholPsi.reset();
    } else {
      cholPsi.reset();
    }

    if (latentDim > 0) {
      const arma::mat eyeLatent = arma::eye(latentDim, latentDim);
      B = arma::solve(IeMinusGamma, eyeLatent, arma::solve_opts::fast);
    } else {
      B.reset();
    }

    Bt = B.t();
    LambdaT = Lambda.t();
    if (alpha.n_cols > 1)
      Rcpp::stop("Alpha matrix with more than one column is not supported.");
    if (alpha.n_elem > 0)
      alphaRow = alpha.t();
    else
      alphaRow = arma::rowvec(latentDim, arma::fill::zeros);

    if (alphaRow.n_elem != latentDim)
      Rcpp::stop("Alpha has incompatible dimensions.");

    const arma::uword manifestDim = Theta.n_rows;
    thetaStd.set_size(manifestDim);
    thetaInvVar.set_size(manifestDim);
    thetaLogNorm.set_size(manifestDim);

    const double log2pi = std::log(2.0 * M_PI);
    for (arma::uword j = 0; j < manifestDim; ++j) {
      const double variance = Theta(j, j);
      const double sd = std::sqrt(variance);
      thetaStd[j] = sd;
      if (sd > 0) {
        const double invVar = 1.0 / (sd * sd);
        thetaInvVar[j] = invVar;
        thetaLogNorm[j] = (-0.5) * log2pi - std::log(sd);
      } else {
        thetaInvVar[j] = arma::datum::nan;
        thetaLogNorm[j] = arma::datum::nan;
      }
    }

    thresholdRows.assign(manifestDim, arma::vec());
    for (arma::uword j = 0; j < manifestDim; ++j) {
      if (j < isordered.n_elem && isordered[j] > 0) {
        const arma::uword idx = isordered[j] - 1L;
        if (idx < Thresholds.n_rows)
          thresholdRows[j] = Thresholds.row(idx).t();
      }
    }

    predictorCache.clear();
    predictorCacheInitialized = false;
  }

  void initializeFreeParams(const Rcpp::List& modFilled) {
    freeLambda.clear();
    freeTau.clear();
    freeGamma.clear();
    freeAlpha.clear();
    freePsi.clear();
    freeTheta.clear();
    freeThresholds.clear();
    freeLabelOrder.clear();

    if (!modFilled.containsElementNamed("free"))
      return;

    Rcpp::List freeR = modFilled["free"];
    std::unordered_set<std::string> seen;

    auto parseMatrix = [&](const char* name, std::vector<FreeEntry>& out) {
      if (!freeR.containsElementNamed(name))
        return;

      Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(freeR[name]);
      if (df.nrows() == 0)
        return;

      std::vector<int> rows = Rcpp::as<std::vector<int>>(df["row"]);
      std::vector<int> cols = Rcpp::as<std::vector<int>>(df["col"]);
      std::vector<int> sym = Rcpp::as<std::vector<int>>(df["symmetric"]);
      std::vector<std::string> labels = Rcpp::as<std::vector<std::string>>(df["label"]);

      const std::size_t len = rows.size();
      out.reserve(len);

      for (std::size_t i = 0; i < len; ++i) {
        const int r = rows[i];
        const int c = cols[i];
        if (r == NA_INTEGER || c == NA_INTEGER)
          continue;
        if (i >= labels.size())
          continue;
        const std::string label = labels[i];
        if (label.empty() || label == "NA")
          continue;

        FreeEntry entry;
        entry.row = r > 0 ? static_cast<arma::uword>(r - 1) : 0;
        entry.col = c > 0 ? static_cast<arma::uword>(c - 1) : 0;
        entry.symmetric = (i < sym.size() && sym[i]);
        entry.label = label;

        out.push_back(entry);

        if (!seen.count(entry.label)) {
          seen.insert(entry.label);
          freeLabelOrder.push_back(entry.label);
        }
      }
    };

    parseMatrix("lambda", freeLambda);
    parseMatrix("tau", freeTau);
    parseMatrix("gamma", freeGamma);
    parseMatrix("alpha", freeAlpha);
    parseMatrix("psi", freePsi);
    parseMatrix("theta", freeTheta);
    parseMatrix("thresholds", freeThresholds);
  }

  void computePredictors(const arma::mat& Zq,
                         arma::mat& latent,
                         arma::mat& S,
                         arma::mat& V) const {
    if (!cholPsiValid)
      Rcpp::stop("Psi is not positive definite; cannot compute predictors.");

    const arma::uword latentDim = B.n_rows;

    if (latentDim == 0) {
      latent.reset();
      S.set_size(N, 0);
      V.set_size(N, k);
      V.each_row() = tau.t();
      return;
    }

    if (cholPsi.n_elem > 0 && Zq.n_cols > 0)
      latent = Zq * cholPsi;
    else
      latent = arma::mat(Zq.n_rows, latentDim, arma::fill::zeros);
    latent.each_row() += alphaRow;

    if (latent.n_rows != N)
      Rcpp::stop("Quadrature nodes have inconsistent number of rows.");

    S = latent * Bt;

    V.set_size(N, k);
    V.each_row() = tau.t();
    if (LambdaT.n_elem > 0)
      V += S * LambdaT;
  }

  void ensurePredictorCacheAllocated() const {
    if (predictorCacheInitialized)
      return;
#ifdef _OPENMP
#pragma omp critical(modsem_gsem_predictor_cache_init)
#endif
    {
      if (!predictorCacheInitialized) {
        predictorCache.assign(Z.size(), PredictorCacheEntry());
        predictorCacheInitialized = true;
      }
    }
  }

  const PredictorCacheEntry& getPredictorCacheEntry(const std::size_t q) const {
    ensurePredictorCacheAllocated();
    if (q >= predictorCache.size())
      Rcpp::stop("Quadrature index out of range.");

    PredictorCacheEntry& entry = predictorCache[q];
    if (!entry.valid) {
#ifdef _OPENMP
#pragma omp critical(modsem_gsem_predictor_cache_fill)
#endif
      {
        if (!entry.valid) {
          arma::mat latent;
          arma::mat S;
          arma::mat V;
          computePredictors(Z[q], latent, S, V);
          entry.latent = std::move(latent);
          entry.S = std::move(S);
          entry.V = std::move(V);
          entry.valid = true;
        }
      }
    }

    return entry;
  }

  Rcpp::NumericVector gradientQ(const arma::mat &P, const int ncores = 1) const {
    if (P.n_rows != N)
      Rcpp::stop("`P` has incompatible number of rows.");
    if (P.n_cols != Z.size())
      Rcpp::stop("`P` has incompatible number of columns.");
    if (!cholPsiValid)
      Rcpp::stop("Psi is not positive definite; cannot compute gradients.");

    std::vector<std::string> labelOrder = freeLabelOrder;
    std::unordered_map<std::string, double> gradMap;
    gradMap.reserve(labelOrder.size());
    for (const auto& label : labelOrder)
      gradMap[label] = 0.0;

    auto addValue = [&](const std::string& label, double value) {
      auto it = gradMap.find(label);
      if (it == gradMap.end()) {
        gradMap.emplace(label, value);
        labelOrder.push_back(label);
      } else {
        it->second += value;
      }
    };

    const int nThreads = resolve_ncores(ncores);
    ensurePredictorCacheAllocated();

    arma::vec gradTau = arma::zeros<arma::vec>(tau.n_rows);
    arma::mat gradLambda = arma::zeros<arma::mat>(Lambda.n_rows, Lambda.n_cols);
    arma::mat gradTheta = arma::zeros<arma::mat>(Theta.n_rows, Theta.n_cols);
    arma::mat gradThresholds = arma::zeros<arma::mat>(Thresholds.n_rows, Thresholds.n_cols);
    arma::mat gradB = arma::zeros<arma::mat>(B.n_rows, B.n_cols);
    arma::mat gradChol = arma::zeros<arma::mat>(cholPsi.n_rows, cholPsi.n_cols);
    arma::vec gradAlphaVec = arma::zeros<arma::vec>(alphaRow.n_elem);

    std::vector<arma::vec> gradTauList(nThreads);
    std::vector<arma::mat> gradLambdaList(nThreads);
    std::vector<arma::mat> gradThetaList(nThreads);
    std::vector<arma::mat> gradThresholdsList(nThreads);
    std::vector<arma::mat> gradBList(nThreads);
    std::vector<arma::mat> gradCholList(nThreads);
    std::vector<arma::vec> gradAlphaList(nThreads);

    for (int t = 0; t < nThreads; ++t) {
      gradTauList[t].zeros(gradTau.n_elem);
      gradLambdaList[t].zeros(gradLambda.n_rows, gradLambda.n_cols);
      gradThetaList[t].zeros(gradTheta.n_rows, gradTheta.n_cols);
      gradThresholdsList[t].zeros(gradThresholds.n_rows, gradThresholds.n_cols);
      gradBList[t].zeros(gradB.n_rows, gradB.n_cols);
      gradCholList[t].zeros(gradChol.n_rows, gradChol.n_cols);
      gradAlphaList[t].zeros(gradAlphaVec.n_elem);
    }

    const double invSqrt2Pi = 1.0 / std::sqrt(2.0 * M_PI);
    const std::size_t nQ = Z.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for (std::size_t q = 0; q < nQ; ++q) {
      int threadId = 0;
#ifdef _OPENMP
      threadId = omp_get_thread_num();
#endif
      arma::vec& gradTauLocal = gradTauList[threadId];
      arma::mat& gradLambdaLocal = gradLambdaList[threadId];
      arma::mat& gradThetaLocal = gradThetaList[threadId];
      arma::mat& gradThresholdsLocal = gradThresholdsList[threadId];
      arma::mat& gradBLocal = gradBList[threadId];
      arma::mat& gradCholLocal = gradCholList[threadId];
      arma::vec& gradAlphaLocal = gradAlphaList[threadId];

      const PredictorCacheEntry& predictors = getPredictorCacheEntry(q);
      const arma::mat& latent = predictors.latent;
      const arma::mat& S = predictors.S;
      const arma::mat& V = predictors.V;
      const auto weights = P.col(q);

      unsigned offset = 0;
      for (int pi = 0; pi < p; ++pi) {
        const int np = n[pi];
        const int end = offset + np - 1;

        const arma::mat& Yp = Y[pi];
        const arma::uvec& colidxp = colIdxPatterns[pi];
        auto weightSub = weights.subvec(offset, end);
        const auto Sblock = S.rows(offset, end);
        arma::mat gradSblock(np, B.n_cols, arma::fill::zeros);
        const arma::uword npU = static_cast<arma::uword>(np);

        arma::vec scoreMean(npU);
        arma::vec scoreVar(npU);
        arma::vec lower(npU);
        arma::vec upper(npU);
        arma::vec phiLower(npU);
        arma::vec phiUpper(npU);
        arma::vec lowerDiff(npU);
        arma::vec upperDiff(npU);
        arma::uvec catIdx(npU);
        arma::uvec catIdxUpper(npU);

        for (arma::uword idx = 0; idx < colidxp.n_elem; ++idx) {
          const arma::uword j = colidxp[idx];
          scoreMean.zeros();
          scoreVar.zeros();

          const arma::vec vj = V.col(j).subvec(offset, end);
          const auto yj = Yp.col(j);

          if (isordered[j]) {
            const arma::vec& thresholdRow = thresholdRows[j];
            if (thresholdRow.n_elem == 0)
              continue;

            const arma::uword threshRowIdx = isordered[j] - 1;
            const double sd = thetaStd[j];
            if (!std::isfinite(sd) || sd <= 0.0)
              continue;

            const arma::uword nThresholds = thresholdRow.n_elem;
            if (nThresholds == 0)
              continue;

            for (arma::uword obs = 0; obs < npU; ++obs) {
              const arma::sword raw = static_cast<arma::sword>(yj[obs]);
              const arma::uword cat = raw > 0 ? static_cast<arma::uword>(raw - 1) : 0;
              catIdx[obs] = cat;
              catIdxUpper[obs] = cat + 1;
              const arma::uword li = std::min(cat, nThresholds - 1);
              const arma::uword ui = std::min(catIdxUpper[obs], nThresholds - 1);
              lower[obs] = thresholdRow[li];
              upper[obs] = thresholdRow[ui];
            }

            arma::vec prob = pnormOrderedProbit(vj, lower, upper, 0, sd, false);

            arma::vec zLower = (lower - vj) / sd;
            arma::vec zUpper = (upper - vj) / sd;

            phiLower = arma::exp(-0.5 * arma::square(zLower)) * invSqrt2Pi;
            phiUpper = arma::exp(-0.5 * arma::square(zUpper)) * invSqrt2Pi;
            lowerDiff = lower - vj;
            upperDiff = upper - vj;

            for (arma::uword obs = 0; obs < npU; ++obs) {
              if (std::isfinite(lower[obs]))
                lowerDiff[obs] = lower[obs] - vj[obs];
              else
                phiLower[obs] = 0.0;

              if (std::isfinite(upper[obs]))
                upperDiff[obs] = upper[obs] - vj[obs];
              else
                phiUpper[obs] = 0.0;
            }

            arma::uvec valid = arma::find(prob > 0);
            if (valid.n_elem == 0)
              continue;

            const double sdCubed = sd * sd * sd;
            scoreMean.elem(valid) = (phiLower.elem(valid) - phiUpper.elem(valid)) /
              (sd * prob.elem(valid));

            scoreVar.elem(valid) = (phiLower.elem(valid) % lowerDiff.elem(valid) -
              phiUpper.elem(valid) % upperDiff.elem(valid)) /
              (2.0 * sdCubed * prob.elem(valid));

            for (arma::uword obs = 0; obs < npU; ++obs) {
              const double w = weightSub[obs];
              if (w == 0.0 || !std::isfinite(prob[obs]) || prob[obs] <= 0.0)
                continue;

              const double invSdProb = 1.0 / (sd * prob[obs]);

              if (std::isfinite(lower[obs]) && catIdx[obs] < gradThresholdsLocal.n_cols)
                gradThresholdsLocal(threshRowIdx, catIdx[obs]) += -phiLower[obs] * invSdProb * w;

              if (std::isfinite(upper[obs]) && catIdxUpper[obs] < gradThresholdsLocal.n_cols)
                gradThresholdsLocal(threshRowIdx, catIdxUpper[obs]) += phiUpper[obs] * invSdProb * w;
            }
          } else {
            const double invVar = thetaInvVar[j];
            if (!std::isfinite(invVar))
              continue;

            const arma::vec residual = yj - vj;
            scoreMean = residual * invVar;

            const double invVarSq = invVar * invVar;
            scoreVar = 0.5 * (arma::square(residual) * invVarSq - invVar);
          }

          arma::vec weightedScore = weightSub % scoreMean;
          const double sumScore = arma::accu(weightedScore);

          if (gradTauLocal.n_elem > j)
            gradTauLocal[j] += sumScore;

          if (gradLambdaLocal.n_cols && Sblock.n_cols) {
            const arma::rowvec ds_dS = Lambda.row(j);
            gradLambdaLocal.row(j) += weightedScore.t() * Sblock;
            gradSblock += weightedScore * ds_dS;
          }

          if (gradThetaLocal.n_rows > j && gradThetaLocal.n_cols > j)
            gradThetaLocal(j, j) += arma::dot(weightSub, scoreVar);
        }

        if (B.n_cols > 0) {
          const arma::mat latentBlock = latent.rows(offset, end);
          gradBLocal += latentBlock.t() * gradSblock;

          const arma::mat gradLatentBlock = gradSblock * B;
          if (gradAlphaLocal.n_elem)
            gradAlphaLocal += arma::sum(gradLatentBlock, 0).t();
          if (cholPsi.n_elem > 0) {
            arma::mat Zblock = Z[q].rows(offset, end);
            if (Zblock.n_cols > cholPsi.n_cols)
              Zblock = Zblock.cols(Zblock.n_cols - cholPsi.n_cols, Zblock.n_cols - 1);
            gradCholLocal += Zblock.t() * gradLatentBlock;
          }
        }

        offset += np;
      }
    }

    for (int t = 0; t < nThreads; ++t) {
      gradTau += gradTauList[t];
      gradLambda += gradLambdaList[t];
      gradTheta += gradThetaList[t];
      gradThresholds += gradThresholdsList[t];
      gradB += gradBList[t];
      gradChol += gradCholList[t];
      gradAlphaVec += gradAlphaList[t];
    }

    arma::mat gradAlphaMat(alpha.n_rows, alpha.n_cols, arma::fill::zeros);
    if (gradAlphaVec.n_elem && alpha.n_elem) {
      if (alpha.n_cols == 1)
        gradAlphaMat.col(0) = gradAlphaVec;
      else
        gradAlphaMat.each_col() = gradAlphaVec;
    }

    arma::mat gradGammaMat(Gamma.n_rows, Gamma.n_cols, arma::fill::zeros);
    if (gradB.n_elem)
      gradGammaMat = (B * gradB * B).t();

    arma::mat gradPsiMat(Psi.n_rows, Psi.n_cols, arma::fill::zeros);
    if (cholPsi.n_elem) {
      arma::mat gradCholUpper = arma::trimatu(gradChol);
      arma::mat temp = arma::solve(arma::trimatu(cholPsi), gradCholUpper, arma::solve_opts::fast);
      arma::mat gradPsiCandidate = 0.5 * temp;
      gradPsiMat = 0.5 * (gradPsiCandidate + gradPsiCandidate.t());
    }

    for (const auto& entry : freeTau) {
      if (entry.row < gradTau.n_elem)
        addValue(entry.label, gradTau[entry.row]);
    }

    for (const auto& entry : freeLambda) {
      if (entry.row < gradLambda.n_rows && entry.col < gradLambda.n_cols)
        addValue(entry.label, gradLambda(entry.row, entry.col));
    }

    for (const auto& entry : freeAlpha) {
      if (entry.row < gradAlphaMat.n_rows && entry.col < gradAlphaMat.n_cols)
        addValue(entry.label, gradAlphaMat(entry.row, entry.col));
    }

    for (const auto& entry : freeTheta) {
      double value = 0.0;
      if (entry.row < gradTheta.n_rows && entry.col < gradTheta.n_cols)
        value = gradTheta(entry.row, entry.col);
      if (entry.symmetric && entry.row != entry.col &&
          entry.col < gradTheta.n_rows && entry.row < gradTheta.n_cols)
        value += gradTheta(entry.col, entry.row);
      addValue(entry.label, value);
    }

    for (const auto& entry : freeThresholds) {
      if (entry.row < gradThresholds.n_rows && entry.col < gradThresholds.n_cols)
        addValue(entry.label, gradThresholds(entry.row, entry.col));
    }

    for (const auto& entry : freeGamma) {
      if (entry.row < gradGammaMat.n_rows && entry.col < gradGammaMat.n_cols)
        addValue(entry.label, gradGammaMat(entry.row, entry.col));
    }

    for (const auto& entry : freePsi) {
      double value = 0.0;

      if (entry.row < gradPsiMat.n_rows && entry.col < gradPsiMat.n_cols) {
        const int scale = entry.symmetric && entry.row != entry.col ? 4L : 1L;
        addValue(entry.label, scale * gradPsiMat(entry.row, entry.col));
      }
    }

    Rcpp::NumericVector out(labelOrder.size());
    for (std::size_t i = 0; i < labelOrder.size(); ++i)
      out[i] = gradMap[labelOrder[i]];
    out.names() = Rcpp::wrap(labelOrder);
    return out;
  }

  arma::mat expectedResponse(const arma::mat Zq) const {
    arma::mat V;
    if (!cholPsiValid) {
      V.set_size(N, k);
      V.fill(arma::datum::nan);
      return V;
    }

    arma::mat latent;
    arma::mat S;
    computePredictors(Zq, latent, S, V);
    return V;
  }

  arma::vec getDensityZq(const std::size_t q, const bool log = false) const {
    if (!cholPsiValid) {
      arma::vec invalid(N);
      invalid.fill(arma::datum::nan);
      return invalid;
    }

    arma::vec ldensity = arma::zeros<arma::vec>(N);
    const arma::mat& V = getPredictorCacheEntry(q).V;

    unsigned offset = 0;
    for (int pi = 0; pi < p; pi++) {
      const int np  = n[pi];
      const int end = offset + np - 1L;

      const arma::mat &Yp = Y[pi];
      const arma::uvec &colidxp = colIdxPatterns[pi];
      const arma::uword npU = static_cast<arma::uword>(np);
      arma::vec lower(npU);
      arma::vec upper(npU);

      for (arma::uword idx = 0; idx < colidxp.n_elem; ++idx) {
        const arma::uword j = colidxp[idx];
        const arma::vec vj = V.col(j).subvec(offset, end);
        const auto yj = Yp.col(j);

        arma::vec ldensj;
        if (isordered[j]) {
          const arma::vec &thresholdsj = thresholdRows[j];
          if (thresholdsj.n_elem == 0) continue;
          const arma::uword nThresholds = thresholdsj.n_elem;
          for (arma::uword obs = 0; obs < npU; ++obs) {
            const arma::sword raw = static_cast<arma::sword>(yj[obs]);
            const arma::uword idxVal = raw > 0 ? static_cast<arma::uword>(raw - 1) : 0;
            const arma::uword lowerIdx = std::min(idxVal, nThresholds - 1);
            const arma::uword upperIdx = std::min(idxVal + 1, nThresholds - 1);
            lower[obs] = thresholdsj[lowerIdx];
            upper[obs] = thresholdsj[upperIdx];
          }
          ldensj = pnormOrderedProbit(vj, lower, upper, 0, thetaStd[j], true);
        } else {
          const arma::vec residual = yj - vj;
          const double invVarHalf = 0.5 * thetaInvVar[j];
          ldensj = thetaLogNorm[j] - invVarHalf * arma::square(residual);
        }

        ldensity.subvec(offset, end) += ldensj;
      }

      offset += np;
    }

    return log ? ldensity : arma::exp(ldensity);
  }

  arma::mat Pi(const bool normalized, const int ncores = 1) {
    if (Z.empty())
      return arma::mat();
    const int nThreads = resolve_ncores(ncores);
    ensurePredictorCacheAllocated();
    arma::mat out(Z[0L].n_rows, Z.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for (int q = 0; q < static_cast<int>(Z.size()); q++)
      out.col(q) = W.col(q) % getDensityZq(static_cast<std::size_t>(q), false);

    return normalized ? out.each_col() / arma::sum(out, 1L): out;
  }

  arma::vec Qi(const arma::mat &P, const int ncores = 1) {
    if (Z.empty())
      return arma::vec();
    const int nThreads = resolve_ncores(ncores);
    ensurePredictorCacheAllocated();
    arma::vec density = arma::zeros<arma::vec>(Z[0L].n_rows);
    std::vector<arma::vec> densityChunks(nThreads);
    for (int t = 0; t < nThreads; ++t)
      densityChunks[t].zeros(density.n_elem);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
    for (int q = 0; q < static_cast<int>(Z.size()); q++) {
      int threadId = 0;
#ifdef _OPENMP
      threadId = omp_get_thread_num();
#endif
      densityChunks[threadId] += P.col(q) % getDensityZq(static_cast<std::size_t>(q), true);
    }

    for (int t = 0; t < nThreads; ++t)
      density += densityChunks[t];

    return density;
  }

  double Q(const arma::mat &P, const int ncores = 1) {
    return arma::accu(Qi(P, ncores)); 
  }

  GSEM_ModelGroup threadClone() const {
    GSEM_ModelGroup c = *this;    // shallow for everything (fast)
                           // Deep-copy ONLY what setParams()/lms_param can modify:
    c.Ie     = arma::mat(Ie);
    c.Lambda = arma::mat(Lambda);
    c.Gamma  = arma::mat(Gamma);
    c.alpha  = arma::mat(alpha);
    c.Psi    = arma::mat(Psi);
    c.Theta  = arma::mat(Theta);

    c.initializeCaches();

    return c;
  }
};


// [[Rcpp::export]]
arma::mat P_Step_GSEM_Group(const Rcpp::List &modelR, const bool normalized, const int ncores) {
  GSEM_ModelGroup M(modelR);
  return M.Pi(normalized, ncores);
}


// [[Rcpp::export]]
double Q_GSEM_Group(const Rcpp::List &modelR, const arma::mat &P, const int ncores) {
  GSEM_ModelGroup M(modelR);
  return M.Q(P, ncores);
}


struct GSEM_Model {
  std::vector<GSEM_ModelGroup> groupModels;
  int ngroups;
  int N = 0;

  explicit GSEM_Model(const Rcpp::List& modFilled) {
    const Rcpp::List groupModelsR = modFilled["models"];
    ngroups = groupModelsR.size();
    groupModels = std::vector<GSEM_ModelGroup>(ngroups);

    N = 0;
    for (int g = 0L; g < ngroups; g++) {
      const Rcpp::List &groupModelR = groupModelsR[g];
      groupModels[g] = GSEM_ModelGroup(groupModelR);
      N += groupModels[g].N;
    }
  }

  arma::mat Pi(const bool normalized, const int ncores = 1) {
    const int zk = groupModels[0L].Z.size();
    arma::mat out(N, zk);
   
    int offset = 0L;
    for (int g = 0L; g < ngroups; g++) {
      const int ng = groupModels[g].N;
      out.rows(offset, offset + ng - 1L) = groupModels[g].Pi(normalized, ncores);
      offset += ng;
    }
   
    return out;
  }

  arma::vec Qi(const arma::mat &P, const int ncores = 1) {
    arma::vec density(N);

    int offset = 0L;
    for (int g = 0L; g < ngroups; g++) {
      const int ng  = groupModels[g].N;
      const int end = offset + ng - 1L;

      const arma::mat &Pg = P.rows(offset, end);

      density.subvec(offset, end) = groupModels[g].Qi(Pg, ncores);
      offset += ng;
    }

    return density;
  }

  double Q(const arma::mat &P, const int ncores = 1) {
    return arma::sum(Qi(P, ncores));
  }

  Rcpp::NumericVector gradientQ(const arma::mat &P, const int ncores = 1) const {
    if (P.n_rows != static_cast<arma::uword>(N))
      Rcpp::stop("`P` has incompatible number of rows.");
    if (ngroups == 0)
      return Rcpp::NumericVector(0);

    std::vector<std::string> labelOrder;
    std::unordered_map<std::string, double> gradMap;

    int offset = 0L;
    for (int g = 0L; g < ngroups; ++g) {
      const int ng = groupModels[g].N;
      const int end = offset + ng - 1L;
      const arma::mat Pg = P.rows(offset, end);

      Rcpp::NumericVector gradGroup = groupModels[g].gradientQ(Pg, ncores);
      Rcpp::CharacterVector names = gradGroup.names();

      for (int i = 0; i < gradGroup.size(); ++i) {
        const std::string label = Rcpp::as<std::string>(names[i]);
        const double value = gradGroup[i];

        auto it = gradMap.find(label);
        if (it == gradMap.end()) {
          gradMap.emplace(label, value);
          labelOrder.push_back(label);
        } else {
          it->second += value;
        }
      }

      offset += ng;
    }

    Rcpp::NumericVector out(labelOrder.size());
    for (std::size_t i = 0; i < labelOrder.size(); ++i)
      out[i] = gradMap[labelOrder[i]];
    out.names() = Rcpp::wrap(labelOrder);
    return out;
  }


  GSEM_Model threadClone() const {
    GSEM_Model c = *this;

    for (int g = 0; g < ngroups; g++)
      c.groupModels[g] = groupModels[g].threadClone();

    return c;
  }
};


// [[Rcpp::export]]
arma::mat P_Step_GSEM(const Rcpp::List &modelR, const bool normalized, const int ncores) {
  GSEM_Model M(modelR);
  return M.Pi(normalized, ncores);
}


// [[Rcpp::export]]
double Q_GSEM(const Rcpp::List &modelR, const arma::mat &P, const int ncores) {
  GSEM_Model M(modelR);
  return M.Q(P, ncores);
}


// [[Rcpp::export]]
arma::vec Qi_GSEM(const Rcpp::List &modelR, const arma::mat &P, const int ncores) {
  GSEM_Model M(modelR);
  return M.Qi(P, ncores);
}


// [[Rcpp::export]]
Rcpp::NumericVector Grad_Q_GSEM_Group(const Rcpp::List &modelR, const arma::mat &P, const int ncores) {
  GSEM_ModelGroup M(modelR);
  return M.gradientQ(P, ncores);
}


// [[Rcpp::export]]
Rcpp::NumericVector Grad_Q_GSEM(const Rcpp::List &modelR, const arma::mat &P, const int ncores) {
  GSEM_Model M(modelR);
  return M.gradientQ(P, ncores);
}


static std::vector<std::string> asStdVector(const Rcpp::CharacterVector& x) {
  std::vector<std::string> out;
  out.reserve(x.size());
  for (auto s : x) out.emplace_back(Rcpp::as<std::string>(s));
  return out;
}


// [[Rcpp::export]]
Rcpp::List expand_quadrature_cpp(const Rcpp::List& Z_list,
                                 const arma::mat& cholPsi,
                                 Rcpp::NumericMatrix gamma,
                                 Rcpp::NumericMatrix omega,
                                 const arma::vec& alphaLatent,
                                 const Rcpp::CharacterVector& xiNames,
                                 const Rcpp::CharacterVector& etaNames,
                                 const Rcpp::CharacterVector& latentNames) {

  const std::size_t latentCount = latentNames.size();
  if (latentCount != static_cast<std::size_t>(cholPsi.n_cols) ||
      latentCount != static_cast<std::size_t>(cholPsi.n_rows))
    Rcpp::stop("Dimension mismatch between `cholPsi` and latent names.");

  Rcpp::List gammaDim = gamma.attr("dimnames");
  Rcpp::List omegaDim = omega.attr("dimnames");
  const Rcpp::CharacterVector gammaRowNames = gammaDim[0];
  const Rcpp::CharacterVector gammaColNames = gammaDim[1];
  const Rcpp::CharacterVector omegaRowNames = omegaDim[0];
  const Rcpp::CharacterVector omegaColNames = omegaDim[1];

  const std::size_t totalVars = gammaColNames.size();
  const std::size_t nXwith = omegaRowNames.size();

  if (totalVars != nXwith + latentCount)
    Rcpp::stop("Unexpected gamma dimension; expected xwith + latent columns.");

  std::vector<std::string> gammaCols = asStdVector(gammaColNames);
  std::vector<std::string> gammaRows = asStdVector(gammaRowNames);
  std::vector<std::string> latentVec = asStdVector(latentNames);
  std::vector<std::string> xwithVec  = asStdVector(omegaRowNames);
  std::vector<std::string> omegaCols = asStdVector(omegaColNames);
  std::vector<std::string> xiVec     = asStdVector(xiNames);
  std::vector<std::string> etaVec    = asStdVector(etaNames);

  std::unordered_map<std::string, std::size_t> gammaColIndex;
  for (std::size_t i = 0; i < gammaCols.size(); ++i)
    gammaColIndex[gammaCols[i]] = i;

  std::unordered_map<std::string, std::size_t> gammaRowIndex;
  for (std::size_t i = 0; i < gammaRows.size(); ++i)
    gammaRowIndex[gammaRows[i]] = i;

  std::unordered_map<std::string, std::size_t> latentIndex;
  for (std::size_t i = 0; i < latentVec.size(); ++i)
    latentIndex[latentVec[i]] = i;

  std::vector<bool> isXi(totalVars, false);
  for (const auto& name : xiVec) {
    auto it = gammaColIndex.find(name);
    if (it != gammaColIndex.end())
      isXi[it->second] = true;
  }

  std::vector<std::vector<std::size_t>> etaParents(etaVec.size());
  std::vector<std::vector<double>>      etaCoefs(etaVec.size());
  for (std::size_t k = 0; k < etaVec.size(); ++k) {
    const auto itRow = gammaRowIndex.find(etaVec[k]);
    if (itRow == gammaRowIndex.end())
      Rcpp::stop("Eta '%s' missing in gamma rows.", etaVec[k]);
    const std::size_t rowIdx = itRow->second;
    Rcpp::NumericMatrix::Row row = gamma.row(rowIdx);
    for (std::size_t j = 0; j < totalVars; ++j) {
      const double val = row[j];
      if (val != 0.0) {
        etaParents[k].push_back(j);
        etaCoefs[k].push_back(val);
      }
    }
  }

  std::vector<std::vector<std::size_t>> xwithParents(xwithVec.size());
  for (std::size_t k = 0; k < xwithVec.size(); ++k) {
    arma::rowvec row = omega.row(static_cast<arma::uword>(k));
    for (std::size_t j = 0; j < omegaCols.size(); ++j) {
      if (row[j] != 0.0) {
        const auto it = gammaColIndex.find(omegaCols[j]);
        if (it == gammaColIndex.end())
          Rcpp::stop("OMEGA column '%s' not found in gamma columns.", omegaCols[j]);
        xwithParents[k].push_back(it->second);
      }
    }
  }

  if (alphaLatent.n_elem != static_cast<int>(latentCount))
    Rcpp::stop("alpha vector length mismatch.");

  Rcpp::List result(Z_list.size());

  for (std::size_t idx = 0; idx < Z_list.size(); ++idx) {
    arma::mat Z = Rcpp::as<arma::mat>(Z_list[idx]);
    if (Z.n_cols != cholPsi.n_cols)
      Rcpp::stop("Quadrature node dimension mismatch with psi.");

    const std::size_t n = Z.n_rows;
    arma::mat Zeta = Z * cholPsi;
    for (std::size_t j = 0; j < latentCount; ++j)
      Zeta.col(j) += alphaLatent[j];

    arma::mat values(n, totalVars, arma::fill::zeros);
    for (std::size_t j = 0; j < latentCount; ++j) {
      const auto it = gammaColIndex.find(latentVec[j]);
      if (it == gammaColIndex.end())
        Rcpp::stop("Latent '%s' missing in gamma columns.", latentVec[j]);
      values.col(it->second) = Zeta.col(j);
    }

    std::vector<bool> defined(totalVars, false);
    for (const auto& name : xiVec) {
      const auto it = gammaColIndex.find(name);
      if (it != gammaColIndex.end())
        defined[it->second] = true;
    }

    arma::mat out(n, nXwith + latentCount, arma::fill::zeros);
    for (std::size_t j = 0; j < latentCount; ++j)
      out.col(nXwith + j) = Z.col(j);

    std::vector<bool> xwithDone(xwithVec.size(), false);
    auto resolveXwith = [&]() -> bool {
      bool progress = false;
      for (std::size_t k = 0; k < xwithVec.size(); ++k) {
        if (xwithDone[k]) continue;
        const auto& parents = xwithParents[k];
        bool ready = true;
        for (const auto parent : parents) {
          if (!defined[parent]) { ready = false; break; }
        }
        if (!ready) continue;

        arma::vec prod(n, arma::fill::ones);
        for (const auto parent : parents)
          prod %= values.col(parent);

        const auto itGamma = gammaColIndex.find(xwithVec[k]);
        if (itGamma == gammaColIndex.end())
          Rcpp::stop("XWITH '%s' missing in gamma columns.", xwithVec[k]);
        const std::size_t tgt = itGamma->second;

        values.col(tgt) = prod;
        defined[tgt] = true;
        out.col(k) = prod;
        xwithDone[k] = true;
        progress = true;
      }
      return progress;
    };

    while (resolveXwith()) {}

    for (std::size_t k = 0; k < etaVec.size(); ++k) {
      const auto itGamma = gammaColIndex.find(etaVec[k]);
      const auto itLatent = latentIndex.find(etaVec[k]);
      if (itGamma == gammaColIndex.end() || itLatent == latentIndex.end())
        Rcpp::stop("Eta '%s' missing in gamma or latent names.", etaVec[k]);
      const std::size_t tgt = itGamma->second;
      const std::size_t latentPos = itLatent->second;

      arma::vec vals = Zeta.col(latentPos);
      const auto& parents = etaParents[k];
      const auto& coefs   = etaCoefs[k];
      for (std::size_t j = 0; j < parents.size(); ++j)
        vals += values.col(parents[j]) * coefs[j];

      values.col(tgt) = vals;
      defined[tgt] = true;

      while (resolveXwith()) {}
    }

    for (std::size_t k = 0; k < xwithVec.size(); ++k) {
      if (!xwithDone[k])
        Rcpp::stop("Failed to resolve latent variable order for quadrature expansion.");
    }

    Rcpp::NumericMatrix outR(n, nXwith + latentCount);
    std::copy(out.begin(), out.end(), outR.begin());

    Rcpp::CharacterVector colNames(nXwith + latentCount);
    for (std::size_t j = 0; j < xwithVec.size(); ++j) colNames[j] = xwithVec[j];
    for (std::size_t j = 0; j < latentVec.size(); ++j) colNames[nXwith + j] = latentVec[j];
    colnames(outR) = colNames;

    result[idx] = outR;
  }

  return result;
}
