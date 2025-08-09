#ifndef mvnorm_h
#define mvnorm_h 

#define BOOST_DISABLE_ASSERTS true

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]

struct ThreadSetter {
  // Saved state
  int         old_threads_mp = 1;
  int         old_threads_blas = 0;       // 0 -> unknown/not saved
  std::string old_vecLib_env;             // for Accelerate
  bool        restore_openblas = false;
  bool        restore_mkl      = false;
  bool        restore_veclib   = false;

  explicit ThreadSetter(const int ncores, const int ncores_blas = 1) {
    // OpenMP
    #ifdef _OPENMP
      if (ncores <= 0)
        Rcpp::stop("ncores must be positive");

      old_threads_mp = omp_get_max_threads();
      omp_set_num_threads(ncores);
    #else
      old_threads_mp = 1;
    #endif

    // BLAS
    #if defined(USE_OPENBLAS)
    
      // OpenBLAS
      extern "C" {
        void openblas_set_num_threads(int);
        int  openblas_get_num_threads(void);
      }
      old_threads_blas = openblas_get_num_threads();
      openblas_set_num_threads(ncores_blas);
      restore_openblas = true;

    #elif defined(USE_MKL)
      // Intel MKL
      extern "C" {
        void mkl_set_num_threads_local(int);
        int  mkl_get_max_threads(void);
      }
      old_threads_blas = mkl_get_max_threads();
      mkl_set_num_threads_local(ncores_blas);
      restore_mkl = true;

    #elif defined(__APPLE__)
      // Apple Accelerate / vecLib: controlled via env var
      if (const char* oldv = std::getenv("VECLIB_MAXIMUM_THREADS")) {
        old_vecLib_env = oldv;  // save to restore later
      } else {
        old_vecLib_env.clear(); // means "unset" originally
      }
      ::setenv("VECLIB_MAXIMUM_THREADS",
          std::to_string(ncores_blas).c_str(), 1);
      restore_veclib = true;

  #else
      // Unknown BLAS; do nothing. Consider documenting how to set its threads.
      (void)ncores_blas;
  #endif
  
  }

  ~ThreadSetter() {
    // Restore BLAS first (in case OpenMP thread count influences it)
    #if defined(USE_OPENBLAS)
      if (restore_openblas && old_threads_blas > 0) {
        extern "C" void openblas_set_num_threads(int);
        openblas_set_num_threads(old_threads_blas);
      }
    #elif defined(USE_MKL)
      if (restore_mkl && old_threads_blas > 0) {
        extern "C" void mkl_set_num_threads_local(int);
        mkl_set_num_threads_local(old_threads_blas);
      }
    #elif defined(__APPLE__)
      if (restore_veclib) {
        if (old_vecLib_env.empty()) {
          // Was unset originally; unset again
          ::unsetenv("VECLIB_MAXIMUM_THREADS");
        } else {
          ::setenv("VECLIB_MAXIMUM_THREADS", old_vecLib_env.c_str(), 1);
        }
      }
    #endif

    // Restore OpenMP
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
