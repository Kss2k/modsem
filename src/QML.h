#ifndef QML_H
#define QML_H


arma::vec muQmlCpp(Rcpp::List m, int t);
arma::vec sigmaQmlCpp(Rcpp::List m, int t);
double varZCpp(arma::mat Omega, arma::mat Sigma1);


#endif 
