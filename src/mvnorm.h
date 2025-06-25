#ifndef mvnorm_h
#define mvnorm_h 

#define BOOST_DISABLE_ASSERTS true

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);

arma::vec dmvnrm_arma_mc(arma::mat const &x,  arma::rowvec const &mean,  
    arma::mat const &sigma, bool const logd);

double totalDmvnWeightedCpp(const arma::vec& mu, const arma::mat& sigma,
    const arma::vec& nu, const arma::mat& S, double tgamma, int n, int d);


arma::vec mahaInt(arma::mat & X, arma::vec & mu, arma::mat & sigma, unsigned int ncores, bool isChol);
arma::vec dmvnInt( arma::mat X, arma::vec mu, arma::mat cholDec, bool log, unsigned int ncores);
arma::vec dmvnfast(arma::mat X,  arma::vec mu,  arma::mat sigma, const bool log, const int ncores, const bool isChol);


#endif // !mvnorm_h
