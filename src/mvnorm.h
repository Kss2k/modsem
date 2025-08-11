#ifndef mvnorm_h
#define mvnorm_h 

#define BOOST_DISABLE_ASSERTS true

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]

struct ThreadSetter {
  int old_threads_mp = 1;

  explicit ThreadSetter(const int ncores) {
    // OpenMP
    #ifdef _OPENMP
      if (ncores <= 0)
        Rcpp::stop("ncores must be positive");

      old_threads_mp = omp_get_max_threads();
      omp_set_num_threads(ncores);
    #else
      old_threads_mp = 1;
    #endif
  }

  ~ThreadSetter() {
    #ifdef _OPENMP
      omp_set_num_threads(old_threads_mp);
    #endif
  }
};


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);

arma::vec dmvnrm_arma_mc(const arma::mat &x, const arma::rowvec &mean,
                         const arma::mat &sigma, const bool log, 
                         const int ncores);

arma::vec rep_dmvnorm(const arma::mat &x, const arma::mat &expected,
                      const arma::mat &sigma, const int t,
                      const int);

double totalDmvnWeightedCpp(const arma::vec& mu, const arma::mat& sigma,
    const arma::vec& nu, const arma::mat& S, double tgamma, int n, int d);


arma::vec mahaInt(arma::mat & X, arma::vec & mu, arma::mat & sigma, unsigned int ncores, bool isChol);
arma::vec dmvnInt( arma::mat X, arma::vec mu, arma::mat cholDec, bool log, unsigned int ncores);
arma::vec dmvnfast(arma::mat X,  arma::vec mu,  arma::mat sigma, const bool log, const int ncores, const bool isChol);


#endif // !mvnorm_h
