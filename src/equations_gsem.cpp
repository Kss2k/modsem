#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>
#include <math.h>

#include "utils.h"
#include "mvnorm.h"
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

  for (arma::uword i = 0; i < n; ++i) {
    const double mu = mean + v[i];
    const double low = lower[i];
    const double upp = upper[i];

    double lower_cdf = 0.0;
    if (!std::isinf(low)) {
      const double z_lower = (low - mu) / (sd * sqrt2);
      lower_cdf = 0.5 * std::erfc(-z_lower);
    }

    double upper_cdf = 1.0;
    if (!std::isinf(upp)) {
      const double z_upper = (upp - mu) / (sd * sqrt2);
      upper_cdf = 0.5 * std::erfc(-z_upper);
    }

    double prob = upper_cdf - lower_cdf;
    if (prob <= 0.0) {
      prob = arma::datum::nan;
    }

    out[i] = log ? std::log(prob) : prob;
  }

  return out;
}


// Gsem model for a single group/cluster
struct GSEM_ModelGroup {
  arma::mat Lambda, tau, Thresholds, Ie, Gamma, Psi, alpha, Theta, W;

  std::vector<arma::mat> Z, Y;
  std::vector<arma::uvec> colIdxPatterns;
  
  arma::uvec n;

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
  }
  
  explicit GSEM_ModelGroup() { // Empty initilizer
    k = 0, N = 0, p = 0;
  }

  arma::mat expectedResponse(const arma::mat Zq) const {
    arma::mat out(N, k);
    arma::vec ivec = arma::ones<arma::vec>(N);

    arma::mat L;
    if (!arma::chol(L, Psi, "lower")) {
      out.fill(arma::datum::nan);
      return out;
    }

    const arma::mat B = arma::inv(Ie - Gamma);
    return ivec * tau.t() + ((ivec * alpha.t() + Zq * L) * B.t()) * Lambda.t();
  }

  arma::vec getDensityZq(const arma::mat Zq, const bool log = false) const {
    arma::vec ldensity = arma::zeros<arma::vec>(N);
    const arma::mat V = expectedResponse(Zq);

    unsigned offset = 0;
    for (int pi = 0; pi < p; pi++) {
      const int np  = n[pi];
      const int end = offset + np - 1L;

      const arma::mat Yp = Y[pi];
      const arma::uvec colidxp = colIdxPatterns[pi];

      for (int idx = 0; idx < colidxp.n_elem; idx++) {
        const unsigned j = colidxp[idx];
        const double sd = std::sqrt(Theta(j, j));

        const arma::vec yj = Yp.col(j);
        const arma::vec vj =  V.col(j).subvec(offset, end);

        arma::vec ldensj;
        if (isordered[j]) {
          const arma::vec thresholdsj = Thresholds.row(isordered[j] - 1L).t();
          const arma::uvec tj = arma::conv_to<arma::uvec>::from(yj) - 1L;
          const arma::vec lower = thresholdsj(tj);
          const arma::vec upper = thresholdsj(tj + 1L);

          if (j == 5) {
          Rcpp::Rcout << "Thresholds:\n" << Thresholds << "\n";
          Rcpp::Rcout << "DATA:\n" << Yp.head_rows(25) << "\n";
          Rcpp::Rcout << "j: " << j << "\n";
          Rcpp::Rcout << "isordered[j]: " << isordered[j] << "\n";
          Rcpp::Rcout << "thresholdsj:\n" << thresholdsj.subvec(0, 25) << "\n";
          Rcpp::Rcout << "tj:\n" << tj.subvec(0, 25) << "\n";
          Rcpp::Rcout << "lower:\n" << lower.subvec(0, 25) << "\n";
          Rcpp::Rcout << "upper:\n" << upper.subvec(0, 25) << "\n";
          }

          ldensj = pnormOrderedProbit(vj, lower, upper, 0, sd);

        } else {
          ldensj = dnorm(yj - vj, 0, sd, true);
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
