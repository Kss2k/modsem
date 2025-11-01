#include <RcppArmadillo.h>
#include <float.h>
#include <cmath>
#include <math.h>

#include "utils.h"
#include "mvnorm.h"
// [[Rcpp::depends(RcppArmadillo)]]


struct GSEMModel {
  arma::mat lambda, tau, thresholds, Ie, gamma, psi, alpha, Y;
  unsigned  k       = 0;
  arma::vec isordered;

  explicit LMSModel(const Rcpp::List& modFilled) {

    Rcpp::List matrices = modFilled["matrices"];
    Rcpp::List info     = modFilled["info"];
    Rcpp::List quad     = modFilled["quad"];

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
  }

  arma::vec expectedResponse(const arma::mat Zi) const {
    arma::vec out(Zi.n_rows);
    arma::vec ivec = arma::ones<arma::vec>(Zi.n_rows);

    const arma::mat L;

    if (!arma::chol(L, Psi, "lower")) {
      out.fill(arma::datum::nan);
      return out;
    }

    const arma::mat B = arma::inv(Ie - gamma);
    return ivec * tau.t() + B * (alpha + Zi * L) * lambda;
  }

  arma::vec getDensityZi(const arma::mat Zi) const {
    arma::vec density = arma::ones<arma::vec>(Zi.n_rows);
    const arma::vec v = expectedResponse(Zi);

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
    for (int q = 0; q < quad.n_elems; q++) 

  }

  LMSModel threadClone() const {
    LMSModel c = *this;    // shallow for everything (fast)
                           // Deep-copy ONLY what setParams()/lms_param can modify:
    c.A     = arma::mat(A);
    c.Oxx   = arma::mat(Oxx);
    c.Oex   = arma::mat(Oex);
    c.Ie    = arma::mat(Ie);
    c.lY    = arma::mat(lY);
    c.lX    = arma::mat(lX);
    c.tY    = arma::mat(tY);
    c.tX    = arma::mat(tX);
    c.Gx    = arma::mat(Gx);
    c.Ge    = arma::mat(Ge);
    c.a     = arma::mat(a);
    c.beta0 = arma::mat(beta0);
    c.Psi   = arma::mat(Psi);
    c.d     = arma::mat(d);
    c.e     = arma::mat(e);

    return c;
  }
};
