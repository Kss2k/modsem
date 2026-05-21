#ifndef mvnorm_h
#define mvnorm_h

#define BOOST_DISABLE_ASSERTS true

#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]

#include "thread_setter.h"  // ThreadSetter (arma-free)


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);

arma::vec dmvnrmArmaMc(const arma::mat &x, const arma::rowvec &mean,
                       const arma::mat &sigma, const bool log, 
                       const int ncores);

arma::vec repDmvnorm(const arma::mat &x, const arma::mat &expected,
                     const arma::mat &sigma, const int t,
                     const int);

double totalDmvnWeighted(const arma::vec& mu,
                         const arma::mat& sigma,
                         const arma::vec& nu,
                         const arma::mat& S,
                         const double tgamma,
                         const int d);

arma::vec mahaInt(arma::mat & X, arma::vec & mu, arma::mat & sigma, unsigned int ncores, bool isChol);
arma::vec dmvnInt( arma::mat X, arma::vec mu, arma::mat cholDec, bool log, unsigned int ncores);
arma::vec dmvnfast(arma::mat X,  arma::vec mu,  arma::mat sigma, const bool log, const int ncores, const bool isChol);


#endif // !mvnorm_h
