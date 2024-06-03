// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <omp.h>
#include "mvnorm.h"


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
                      const arma::mat &sigma, const int t,
                      const int cores) {
  arma::vec out = arma::vec(t);
  int ncolSigma = sigma.n_cols;
  int firstRow, lastRow, firstCol, lastCol;
  for (int i = 0; i < t; i++) {
    firstRow = i * ncolSigma;
    lastRow = (i + 1) * ncolSigma - 1;
    firstCol = 0;
    lastCol = ncolSigma - 1;

    out(i) = dmvnrm_arma_mc(x.row(i), expected.row(i), 
        sigma.submat(firstRow, firstCol, lastRow, lastCol), true, cores)(0);
  }
  return out;
}


// [[Rcpp::export]]
arma::vec dmvnrm_arma_mc(arma::mat const &x,  
                         arma::rowvec const &mean,  
                         arma::mat const &sigma, 
                         bool const logd = true,
                         int const cores = 2) {  
    using arma::uword;
    omp_set_num_threads(cores);
    uword const n = x.n_rows, 
             xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
                constants = -(double)xdim/2.0 * log2pi, 
              other_terms = rootisum + constants;
    
    arma::rowvec z;
    #pragma omp parallel for schedule(static) private(z)
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);   
        out(i) = other_terms - 0.5 * arma::dot(z, z);     
    }  
      
    if (logd)
      return out;
    return exp(out);
}
