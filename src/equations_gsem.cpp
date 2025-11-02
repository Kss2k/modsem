#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>
#include <math.h>

#include "utils.h"
#include "mvnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]


// Gsem model for a single group/cluster
struct GSEMModelGroup {
  arma::mat Lambda, tau, Thresholds, Ie, Gamma, Psi, alpha, Y, Theta;

  std::vector<arma::mat> Z, Y, P;
  std::vector<arma::uvec> colIdxPatterns;

  unsigned k = 0, n = 0, p = 0;
  arma::vec isordered;

  explicit GSEMModelGroup(const Rcpp::List& modFilled) {

    Rcpp::List matrices  = modFilled["matrices"];
    Rcpp::List info      = modFilled["info"];

    Rcpp::List quad      = modFilled["quad"];
    Rcpp::List quad_i    = quad["quad_i"];

    Rcpp::List dataR     = modFilled["data"];
    Rcpp::List dataSplit = dataR["data.split"];
    Rcpp::List colidxR   = dataR["colidx"];

    n = Rcpp::as<unsigned>(data["n"]); // Rows
    k = Rcpp::as<unsigned>(data["k"]); // Cols
    p = Rcpp::as<unsigned>(data["npatterns"]); // patterns

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
    P          = Rcpp::as<arma::mat)(quad["P"]);

    Z = std::vector(Rcpp::length(quad_i));
    for (int q = 0; q < Rcpp::length(quad_i); q++)
      Z[q] = Rcpp::as<arma::mat>(quad_i[q]);

    Y              = std::vector(p);
    colInxPatterns = std::vector(p);

    for (int pi = 0; pi < p; pi++) {
      Y[pi]          = Rcpp::as<arma::mat>(dataSplit[pi]);
      colIdxPatterns = Rcpp::as<arma::uvec>(colidxR[pi]);
    }
  }

  arma::vec expectedResponse(const arma::mat Zq) const {
    arma::vec out(Zq.n_rows);
    arma::vec ivec = arma::ones<arma::vec>(Zq.n_rows);

    const arma::mat L;

    if (!arma::chol(L, Psi, "lower")) {
      out.fill(arma::datum::nan);
      return out;
    }

    const arma::mat B = arma::inv(Ie - gamma);
    return ivec * tau.t() + B * (alpha + Zq * L) * Lambda;
  }

  arma::vec getDensityZq(const arma::mat Zq, const bool log = false) const {
    arma::vec density(Zq.n_rows);
    const arma::vec v = expectedResponse(Zq);

    unsigned offset = 0;
    for (int pi = 0; pi < p; pi++) {
      const arma::mat Yp = Y[p];
      const arma::uvec colidxp = colIdxPatterns[pi];
      const unsigned end = offset + Yp.n_rows - 1L;

      for (auto j = colidxp.begin(); j != colidxp.end(); j++) {
        const double sd = std::sqrt(Theta(j, j));

        const arma::vec yj = Yp.col(j).subvec(offset, end);
        const arma::vec vj =  v.col(j).subvec(offset, end);

        density.subvec(offset, end) = Rcpp::dnorm(yj - vj, 0, sd, log);
      }

      offset = end + 1L;
    }
  }

  arma::vec Qi() {
    arma::vec density = arma::zeros<arma::vec>(Z[0L].n_rows);

    for (int q = 0; q < Z.size(); q++)
      density += P.col(q) * getDensityZq(Z[q], true);

    return density;
  }


  double Q(const bool log = true) {
    return arma::sum(Qi()); 
  }

  LMSModel threadClone() const {
    LMSModel c = *this;    // shallow for everything (fast)
                           // Deep-copy ONLY what setParams()/lms_param can modify:
    c.Ie    = arma::mat(Ie);
    c.Lamba = arma::mat(Lambda);
    c.Gamma = arma::mat(Gamma);
    c.alpha = arma::mat(alpha);
    c.Psi   = arma::mat(Psi);
    c.Theta = arma::mat(Theta);

    return c;
  }
};
