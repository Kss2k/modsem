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
// [[Rcpp::depends(RcppArmadillo)]]


arma::vec dnorm(const arma::vec x, const double mu, const double sd, const bool log = false) {
  const arma::vec e = x - mu;
  const arma::vec ldens = (- 0.5) * std::log(2 * M_PI) - std::log(sd) - (e % e)/(2 * sd * sd);
  return log ? ldens: arma::exp(ldens);
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
      cholPsiValid = arma::chol(cholPsi, Psi, "lower");
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

      Rcpp::DataFrame df = freeR[name];
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

  Rcpp::NumericVector gradientQ(const arma::mat &P) const {
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

    arma::vec gradTau = arma::zeros<arma::vec>(tau.n_rows);
    arma::mat gradLambda = arma::zeros<arma::mat>(Lambda.n_rows, Lambda.n_cols);
    arma::mat gradTheta = arma::zeros<arma::mat>(Theta.n_rows, Theta.n_cols);
    arma::mat gradThresholds = arma::zeros<arma::mat>(Thresholds.n_rows, Thresholds.n_cols);
    arma::mat gradB = arma::zeros<arma::mat>(B.n_rows, B.n_cols);
    arma::mat gradChol = arma::zeros<arma::mat>(cholPsi.n_rows, cholPsi.n_cols);
    arma::vec gradAlphaVec = arma::zeros<arma::vec>(alphaRow.n_elem);

    const double invSqrt2Pi = 1.0 / std::sqrt(2.0 * M_PI);
    arma::mat latent;
    arma::mat S;
    arma::mat V;

    for (std::size_t q = 0; q < Z.size(); ++q) {
      computePredictors(Z[q], latent, S, V);
      const arma::vec weights = P.col(q);

      unsigned offset = 0;
      for (int pi = 0; pi < p; ++pi) {
        const int np = n[pi];
        const int end = offset + np - 1;

        const arma::mat& Yp = Y[pi];
        const arma::uvec& colidxp = colIdxPatterns[pi];
        const arma::vec weightSub = weights.subvec(offset, end);
        const auto Sblock = S.rows(offset, end);
        arma::mat gradSblock(np, B.n_cols, arma::fill::zeros);

        for (arma::uword idx = 0; idx < colidxp.n_elem; ++idx) {
          const arma::uword j = colidxp[idx];
          arma::vec scoreMean(np, arma::fill::zeros);
          arma::vec scoreVar(np, arma::fill::zeros);

          const arma::vec vj = V.col(j).subvec(offset, end);
          const arma::vec yj = Yp.col(j);

          if (isordered[j]) {
            const arma::vec& thresholdRow = thresholdRows[j];
            if (thresholdRow.n_elem == 0)
              continue;

            const arma::uword threshRowIdx = isordered[j] - 1;
            const double sd = thetaStd[j];
            if (!std::isfinite(sd) || sd <= 0.0)
              continue;

            arma::uvec catIdx = arma::conv_to<arma::uvec>::from(yj) - 1;
            arma::uvec catIdxUpper = catIdx + 1;

            arma::vec lower(np, arma::fill::zeros);
            arma::vec upper(np, arma::fill::zeros);
            const arma::uword nThresholds = thresholdRow.n_elem;
            if (nThresholds == 0)
              continue;
            for (arma::uword obs = 0; obs < static_cast<arma::uword>(np); ++obs) {
              const arma::uword li = catIdx[obs] < nThresholds ? catIdx[obs] : nThresholds - 1;
              const arma::uword ui = catIdxUpper[obs] < nThresholds ? catIdxUpper[obs] : nThresholds - 1;
              lower[obs] = thresholdRow[li];
              upper[obs] = thresholdRow[ui];
            }

            arma::vec logProb = pnormOrderedProbit(vj, lower, upper, 0, sd, true);
            arma::vec prob = arma::exp(logProb);

            arma::vec zLower = (lower - vj) / sd;
            arma::vec zUpper = (upper - vj) / sd;

            arma::vec phiLower = arma::exp(-0.5 * arma::square(zLower)) * invSqrt2Pi;
            arma::vec phiUpper = arma::exp(-0.5 * arma::square(zUpper)) * invSqrt2Pi;

            arma::vec lowerDiff(np, arma::fill::zeros);
            arma::vec upperDiff(np, arma::fill::zeros);

            for (arma::uword obs = 0; obs < static_cast<arma::uword>(np); ++obs) {
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

            for (arma::uword obs = 0; obs < static_cast<arma::uword>(np); ++obs) {
              const double w = weightSub[obs];
              if (w == 0.0 || !std::isfinite(prob[obs]) || prob[obs] <= 0.0)
                continue;

              const double invSdProb = 1.0 / (sd * prob[obs]);

              if (std::isfinite(lower[obs]) && catIdx[obs] < gradThresholds.n_cols)
                gradThresholds(threshRowIdx, catIdx[obs]) += -phiLower[obs] * invSdProb * w;

              if (std::isfinite(upper[obs]) && catIdxUpper[obs] < gradThresholds.n_cols)
                gradThresholds(threshRowIdx, catIdxUpper[obs]) += phiUpper[obs] * invSdProb * w;
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

          if (gradTau.n_elem > j)
            gradTau[j] += sumScore;

          if (gradLambda.n_cols && Sblock.n_cols) {
            gradLambda.row(j) += weightedScore.t() * Sblock;
            gradSblock += weightedScore * Lambda.row(j);
          }

          if (gradTheta.n_rows > j && gradTheta.n_cols > j)
            gradTheta(j, j) += arma::dot(weightSub, scoreVar);
        }

        if (B.n_cols > 0) {
          const arma::mat latentBlock = latent.rows(offset, end);
          gradB += latentBlock.t() * gradSblock;

          const arma::mat gradLatentBlock = gradSblock * B;
          if (gradAlphaVec.n_elem)
            gradAlphaVec += arma::sum(gradLatentBlock, 0).t();
          if (cholPsi.n_elem > 0)
            gradChol += Z[q].rows(offset, end).t() * gradLatentBlock;
        }

        offset += np;
      }
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
      gradGammaMat = B.t() * gradB * B.t();

    arma::mat gradPsiMat(Psi.n_rows, Psi.n_cols, arma::fill::zeros);
    if (cholPsi.n_elem) {
      arma::mat lowerGrad = arma::trimatl(gradChol);
      arma::mat halfGrad = 0.5 * lowerGrad;
      arma::mat temp = arma::solve(arma::trimatu(cholPsi.t()), halfGrad.t(), arma::solve_opts::fast);
      gradPsiMat = arma::symmatu(temp.t());
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

    for (const auto& entry : freeGamma)
      if (entry.row < gradGammaMat.n_rows && entry.col < gradGammaMat.n_cols)
        addValue(entry.label, gradGammaMat(entry.row, entry.col));
    for (const auto& entry : freePsi) {
      double value = 0.0;
      if (entry.row < gradPsiMat.n_rows && entry.col < gradPsiMat.n_cols)
        value += gradPsiMat(entry.row, entry.col);
      if (entry.symmetric && entry.row != entry.col &&
          entry.col < gradPsiMat.n_rows && entry.row < gradPsiMat.n_cols)
        value += gradPsiMat(entry.col, entry.row);
      addValue(entry.label, value);
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

  arma::vec getDensityZq(const arma::mat Zq, const bool log = false) const {
    arma::vec ldensity = arma::zeros<arma::vec>(N);
    const arma::mat V = expectedResponse(Zq);

    unsigned offset = 0;
    for (int pi = 0; pi < p; pi++) {
      const int np  = n[pi];
      const int end = offset + np - 1L;

      const arma::mat &Yp = Y[pi];
      const arma::uvec &colidxp = colIdxPatterns[pi];

      for (arma::uword idx = 0; idx < colidxp.n_elem; ++idx) {
        const arma::uword j = colidxp[idx];
        const arma::vec vj = V.col(j).subvec(offset, end);
        const arma::vec yj = Yp.col(j);

        arma::vec ldensj;
        if (isordered[j]) {
          const arma::vec &thresholdsj = thresholdRows[j];
          if (thresholdsj.n_elem == 0) continue;
          const arma::uword nThresholds = thresholdsj.n_elem;
          arma::uvec tj = arma::conv_to<arma::uvec>::from(yj);
          arma::vec lower(np, arma::fill::zeros);
          arma::vec upper(np, arma::fill::zeros);
          for (arma::uword obs = 0; obs < static_cast<arma::uword>(np); ++obs) {
            const arma::uword idx = tj[obs] > 0 ? tj[obs] - 1 : 0;
            const arma::uword lowerIdx = std::min(idx, nThresholds - 1);
            const arma::uword upperIdx = std::min(idx + 1, nThresholds - 1);
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

  arma::mat Pi(const bool normalized) {
    arma::mat out(Z[0L].n_rows, Z.size());
     
    for (int q = 0; q < Z.size(); q++)
      out.col(q) = W.col(q) % getDensityZq(Z[q], false);

    return normalized ? out.each_col() / arma::sum(out, 1L): out;
  }

  arma::vec Qi(const arma::mat &P) {
    arma::vec density = arma::zeros<arma::vec>(Z[0L].n_rows);

    for (int q = 0; q < Z.size(); q++)
      density += P.col(q) % getDensityZq(Z[q], true);

    return density;
  }

  double Q(const arma::mat &P) {
    return arma::accu(Qi(P)); 
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
arma::mat P_Step_GSEM_Group(const Rcpp::List &modelR, const bool normalized) {
  GSEM_ModelGroup M(modelR);
  return M.Pi(normalized);
}


// [[Rcpp::export]]
double Q_GSEM_Group(const Rcpp::List &modelR, const arma::mat &P) {
  GSEM_ModelGroup M(modelR);
  return M.Q(P);
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

  arma::mat Pi(const bool normalized) {
    const int zk = groupModels[0L].Z.size();
    arma::mat out(N, zk);
   
    int offset = 0L;
    for (int g = 0L; g < ngroups; g++) {
      const int ng = groupModels[g].N;
      out.rows(offset, offset + ng - 1L) = groupModels[g].Pi(normalized);
      offset += ng;
    }
   
    return out;
  }

  arma::vec Qi(const arma::mat &P) {
    arma::vec density(N);

    int offset = 0L;
    for (int g = 0L; g < ngroups; g++) {
      const int ng  = groupModels[g].N;
      const int end = offset + ng - 1L;

      const arma::mat &Pg = P.rows(offset, end);

      density.subvec(offset, end) = groupModels[g].Qi(Pg);
      offset += ng;
    }

    return density;
  }

  double Q(const arma::mat &P) {
    return arma::sum(Qi(P));
  }

  Rcpp::NumericVector gradientQ(const arma::mat &P) const {
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

      Rcpp::NumericVector gradGroup = groupModels[g].gradientQ(Pg);
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
arma::mat P_Step_GSEM(const Rcpp::List &modelR, const bool normalized) {
  GSEM_Model M(modelR);
  return M.Pi(normalized);
}


// [[Rcpp::export]]
double Q_GSEM(const Rcpp::List &modelR, const arma::mat &P) {
  GSEM_Model M(modelR);
  return M.Q(P);
}


// [[Rcpp::export]]
arma::vec Qi_GSEM(const Rcpp::List &modelR, const arma::mat &P) {
  GSEM_Model M(modelR);
  return M.Qi(P);
}


// [[Rcpp::export]]
Rcpp::NumericVector Grad_Q_GSEM_Group(const Rcpp::List &modelR, const arma::mat &P) {
  GSEM_ModelGroup M(modelR);
  return M.gradientQ(P);
}


// [[Rcpp::export]]
Rcpp::NumericVector Grad_Q_GSEM(const Rcpp::List &modelR, const arma::mat &P) {
  GSEM_Model M(modelR);
  return M.gradientQ(P);
}
