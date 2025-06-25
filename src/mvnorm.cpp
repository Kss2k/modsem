// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "mvnorm.h"
#include <cfloat>


static double const log2pi = std::log(2.0 * M_PI);


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}


// [[Rcpp::export]]
arma::vec rep_dmvnorm(const arma::mat &x, 
                      const arma::mat &expected, 
                      const arma::mat &sigma, const int t) {
  arma::vec out = arma::vec(t);
  int ncolSigma = sigma.n_cols;
  int firstRow, lastRow, firstCol, lastCol;
  for (int i = 0; i < t; i++) {
    firstRow = i * ncolSigma;
    lastRow = (i + 1) * ncolSigma - 1;
    firstCol = 0;
    lastCol = ncolSigma - 1;

    out(i) = dmvnrm_arma_mc(x.row(i), expected.row(i), 
        sigma.submat(firstRow, firstCol, lastRow, lastCol), true)(0);
  }
  return out;
}


// [[Rcpp::export]]
arma::vec dmvnrm_arma_mc(arma::mat const &x,  
                         arma::rowvec const &mean,  
                         arma::mat const &sigma, 
                         bool const logd = true) {
    using arma::uword;
    uword const n = x.n_rows, 
             xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
                constants = -(double)xdim/2.0 * log2pi, 
              other_terms = rootisum + constants;
    
    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);   
        out(i) = other_terms - 0.5 * arma::dot(z, z);     
    }  
      
    if (logd)
      return out;
    return exp(out);
}


// [[Rcpp::export]]
double totalDmvnWeightedCpp(const arma::vec& mu,
                            const arma::mat& sigma,
                            const arma::vec& nu,
                            const arma::mat& S,
                            double tgamma,
                            int n,
                            int d) {
  if (!sigma.is_finite()) return NA_REAL;

  if (arma::any(arma::diagvec(sigma) <= 0)) return NA_REAL;

  arma::mat L;
  if (!arma::chol(L, sigma, "lower")) return NA_REAL;

  double log_det_sigma = 2.0 * arma::sum(log(L.diag()));
  arma::mat sigma_inv = arma::inv_sympd(sigma);
  if (!sigma_inv.is_finite()) return NA_REAL;

  arma::vec diff = nu - mu;
  double mahalanobis_term = tgamma * arma::as_scalar(diff.t() * sigma_inv * diff);

  double trace_term = arma::accu(sigma_inv % S);

  double log_likelihood = -0.5 * (
    tgamma * d * std::log(2.0 * M_PI) +
    tgamma * log_det_sigma +
    trace_term +
    mahalanobis_term
  );

  return log_likelihood;
}


arma::mat dmvnfast(arma::mat X,  
                   arma::mat mu,  
                   arma::mat sigma, 
                   const bool log,
                   const int ncores,
                   const bool isChol) { 
  using namespace Rcpp;
  // Code is copied from mvnfast package
  // I want to call it directly from C++, not from R
  // Copyright (C) 2014 Matteo Fasiolo  matteo.fasiolo@gmail.com
  
  try{
    double df = -1.0;
    
    if(ncores == 0) stop("ncores has to be positive.");
    if (X.n_cols != mu.n_elem) Rcpp::stop("X.n_cols != mu.n_elem"); 
    if (X.n_cols != sigma.n_cols) Rcpp::stop("X.n_cols != sigma.n_cols"); 
    if (sigma.n_rows != sigma.n_cols) Rcpp::stop("sigma.n_rows != sigma.n_cols"); 
    
    // Here we set the number of OMP threads, but before we save the original
    // number of threads, so we can re-set before returning.
    int ncores_0;
    #ifdef _OPENMP
    #pragma omp parallel num_threads(1)
    {
    #pragma omp single
     ncores_0 = omp_get_num_threads();
    }
    omp_set_num_threads(ncores);
    #endif
    
    // Calculate cholesky dec. unless sigma is alread a cholesky dec.
    arma::mat cholDec;
    if (!isChol) cholDec = arma::chol(sigma);
    else cholDec = sigma;
    
    // Dropping the dimensionality of the output vector
    const arma::mat out = dmvtInt( X, mu, cholDec, log, df, ncores);
    
    #ifdef _OPENMP
    omp_set_num_threads(ncores_0);
    #endif
    
    return out;
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);

  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return arma::mat();
}


arma::vec dmvtInt( arma::mat X, arma::vec mu, arma::mat cholDec, bool log, double df, unsigned int ncores) {
  using namespace arma;
  
  unsigned int d = X.n_cols;
  
  vec out = mahaInt(X, mu, cholDec, ncores, true);
  
  if( df <= 0.0 ){ // Multivariate normal density OR ...
    
    out = - 0.5 * out - ( (d / 2.0) * std::log(2.0 * M_PI) + sum(arma::log(cholDec.diag())) );
    
  } else { // ... multivariate Student-t density
    
  #ifdef _OPENMP
  #pragma omp parallel num_threads(ncores) if(ncores > 1)
  {
  #endif
  
  uint32_t ii;  
  uint32_t n = X.n_rows;  
  double logDet = sum(arma::log(cholDec.diag())); 
  double c = lgamma((d + df)/2.0) - (lgamma(df/2.0) + logDet + d/2.0 * std::log(M_PI * df));

  #ifdef _OPENMP
  #pragma omp for schedule(static)
  #endif
  for(ii = 0; ii < n; ii++) {
     out.at(ii) = c - 0.5 * (df + d) * log1p(out.at(ii)/df);
  }
    
  #ifdef _OPENMP
  }
  #endif
    
  }
  
  if (log == false) out = exp(out);
  
  return( out );
}


/* 
 *  Internal C++ function for Mahalanobis distance
*/
arma::vec mahaInt(arma::mat & X,  
                  arma::vec & mu,  
                  arma::mat & sigma,
                  const unsigned int ncores,
                  const bool isChol = false) {
  using namespace arma;
  
  // Some sanity checks 
  if(ncores == 0) Rcpp::stop("ncores has to be positive.");
  if(mu.n_elem != sigma.n_cols) Rcpp::stop("The mean vector has a different dimensions from the covariance matrix.");
  if(X.n_cols != sigma.n_cols)  Rcpp::stop("The number of columns of X is different from the dimension of the covariance matrix.");
                   
  // Calculate transposed cholesky dec. unless sigma is alread a cholesky dec.
  mat cholDec;
  if( isChol == false ) {
     cholDec = trimatl(chol(sigma).t());
  }
  else{
     cholDec = trimatl(sigma.t()); 
     if(any(cholDec.diag() <= 0.0))  Rcpp::stop("The supplied cholesky decomposition has values <= 0.0 on the main diagonal.");
  }
  
  vec D = cholDec.diag();
    
  vec out(X.n_rows);
  
  #ifdef _OPENMP
  #pragma omp parallel num_threads(ncores) if(ncores > 1)                       
  {
  #endif
  
  // Declaring some private variables
  uint32_t d = X.n_cols;
  uint32_t n = X.n_rows;
  
  vec tmp(d);  
    
  double acc;
  uint32_t icol, irow, ii;  
  
  // For each of the "n" random vectors, forwardsolve the corresponding linear system.
  // Forwardsolve because I'm using the lower triangle Cholesky.
  #ifdef _OPENMP
  #pragma omp for schedule(static)
  #endif
  for(icol = 0; icol < n; icol++)
  {
        
    for(irow = 0; irow < d; irow++)
    {
     acc = 0.0;
     
     for(ii = 0; ii < irow; ii++) acc += tmp.at(ii) * cholDec.at(irow, ii);
     
     tmp.at(irow) = ( X.at(icol, irow) - mu.at(irow) - acc ) / D.at(irow);
    }
    
    out.at(icol) = sum(square(tmp)); 
  }
  
  #ifdef _OPENMP
  }
  #endif
    
  return out;
}
