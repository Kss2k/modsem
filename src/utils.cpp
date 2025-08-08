#include <RcppArmadillo.h>
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]


std::vector<arma::vec> as_vec_of_vec(const Rcpp::List& L) {
  const std::size_t J = L.size();
  std::vector<arma::vec> out;
  out.reserve(J);

  for (std::size_t j = 0; j < J; ++j)
    out.emplace_back(Rcpp::as<arma::vec>(L[j]));

  return out;
}


std::vector<std::vector<arma::vec>> as_vec_of_vec_of_vec(const Rcpp::List& L) {
  const std::size_t J = L.size();
  std::vector<std::vector<arma::vec>> out;
  out.reserve(J);

  for (std::size_t j = 0; j < J; ++j)
    out.emplace_back(as_vec_of_vec(L[j]));

  return out;
}


std::vector<arma::mat> as_vec_of_mat(const Rcpp::List& L) {
  const std::size_t J = L.size();
  std::vector<arma::mat> out;
  out.reserve(J);

  for (std::size_t j = 0; j < J; ++j)
    out.emplace_back(Rcpp::as<arma::mat>(L[j]));

  return out;
}


std::vector<std::vector<arma::mat>> as_vec_of_vec_of_mat(const Rcpp::List& L) {
  const std::size_t J = L.size();
  std::vector<std::vector<arma::mat>> out;
  out.reserve(J);

  for (std::size_t j = 0; j < J; ++j)
    out.emplace_back(as_vec_of_mat(L[j]));

  return out;
}


std::vector<arma::uvec> as_vec_of_uvec(const Rcpp::List& L) {
  const std::size_t J = L.size();
  std::vector<arma::uvec> out;
  out.reserve(J);

  for (std::size_t j = 0; j < J; ++j)
    out.emplace_back(Rcpp::as<arma::uvec>(L[j]));

  return out;
}
