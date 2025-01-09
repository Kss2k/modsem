// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>


// [[Rcpp::export]]
arma::vec calcSESimpleSlopes(arma::mat const &X, arma::mat const &V) {
  arma::vec s = arma::vec(X.n_rows);

  for (int i = 0; i < (int)(X.n_rows); i++)
    s(i) = arma::as_scalar(X.row(i) * V * X.row(i).t());

  return s;
}
