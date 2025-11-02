#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>
#include <math.h>

#include "utils.h"
#include "mvnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]


struct GSEMModel {
  arma::mat lambda, tau, thresholds, Ie, gamma, psi, alpha, Y;
  std::vector<arma::mat> Z;
  unsigned  k       = 0;
  arma::vec isordered;

  explicit LMSModel(const Rcpp::List& modFilled) {

    Rcpp::List matrices = modFilled["matrices"];
    Rcpp::List info     = modFilled["info"];
    Rcpp::List quad     = modFilled["quad"];
    Rcpp::List quad_i   = modFilled["quad_i"];

    k  = Rcpp::as<unsigned>(quad["k"]);

    // one-liners, no loops
    Ie         = Rcpp::as<arma::mat>(matrices["Ieta"]);
    Lambda     = Rcpp::as<arma::mat>(matrices["lambda"]);
    Tau        = Rcpp::as<arma::mat>(matrices["tau"]);
    Gamma      = Rcpp::as<arma::mat>(matrices["gamma"]);
    alpha      = Rcpp::as<arma::mat>(matrices["alpha"]);
    Psi        = Rcpp::as<arma::mat>(matrices["psi"]);
    Theta      = Rcpp::as<arma::mat>(matrices["theta"]);
    Thresholds = Rcpp::as<arma::mat>(matrices["thresholds"]);
    Y          = Rcpp::as<arma::mat>(matrices["data"]);
    isordered  = Rcpp::as<arma::vec>(matrices["isordered"]);


    Z = std::vector(Rcpp::length(quad_i));
    for (int q = 0; q < Rcpp::length(quad_i); q++)
      Z[q] = Rcpp::as<arma::mat>(quad_i[q]);
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
    return ivec * tau.t() + B * (alpha + Zq * L) * lambda;
  }

  arma::vec getDensityZq(const arma::mat Zq) const {
    arma::vec density = arma::ones<arma::vec>(Zq.n_rows);
    const arma::vec v = expectedResponse(Zq);

    for (int j = 0; j < Y.n_cols; j++) {
      const double sd = std::sqrt(Theta(j, j));

      for (int i = 0; i < Y.n_rows; i++) {
        const double yij = Y(i, j);

        if (std::isnan(yij))
          continue;

        density(i) = Rcpp::dnorm(yij, vij, sd, true);
      }
    }
  }

  arma::vec Qi() {
    arma::vec density = arma::ones<arma::vec>(Z[q].n_rows);

    for (int q = 0; q < quad.n_elems; q++)
      density *= getDensityZq(Z[q]);

    return density
  }


  arma::vec Q() {

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
