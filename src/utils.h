#ifndef UTILS_H
#define UTILS_H

std::vector<arma::vec> as_vec_of_vec(const Rcpp::List& L);
std::vector<std::vector<arma::vec>> as_vec_of_vec_of_vec(const Rcpp::List& L);
std::vector<arma::mat> as_vec_of_mat(const Rcpp::List& L);
std::vector<std::vector<arma::mat>> as_vec_of_vec_of_mat(const Rcpp::List& L);
std::vector<arma::uvec> as_vec_of_uvec(const Rcpp::List& L);

#endif
