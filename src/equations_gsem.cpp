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


// Gsem model for a single group/cluster
struct GSEM_ModelGroup {
  arma::mat Lambda, tau, Thresholds, Ie, Gamma, Psi, alpha, Theta, W;

  std::vector<arma::mat> Z, Y;
  std::vector<arma::uvec> colIdxPatterns;
  
  arma::uvec n;

  unsigned k = 0, p = 0, N = 0;
  arma::vec isordered;

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
    isordered  = Rcpp::as<arma::vec>(matrices["isordered"]);

    W = Rcpp::as<arma::mat>(quad["w"]);
    Z = std::vector<arma::mat>(quad_i.length());

    for (int q = 0; q < quad_i.length(); q++) {
      Z[q] = Rcpp::as<arma::mat>(quad_i[q]);
    }

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
    arma::vec density(N);
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

        density.subvec(offset, end) = dnorm(yj - vj, 0, sd, log);
      }

      offset += np;
    }

    return density;
  }

  arma::mat Pi(const bool normalized) {
    arma::mat out(Z[0L].n_rows, Z.size());
     
    for (int q = 0; q < Z.size(); q++) {
      Rcpp::Rcout << "W.col(" << g << "):\n" << W.head_rows(5).col(q) << "\n";
      Rcpp::Rcout << "density(Zq[" << q <<"]):\n" << getDensityZq(Z[q], false).subvec(0, 4) << "\n";
      out.col(q) = W.col(q) % getDensityZq(Z[q], false);
    }

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
