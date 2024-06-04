#ifndef QML_H
#define QML_H


arma::mat muQmlCpp(Rcpp::List m, int t);
arma::mat sigmaQmlCpp(Rcpp::List m, int t);
arma::mat varZCpp(arma::mat Omega, arma::mat Sigma1, int numEta);
double varZSubOmega(arma::mat Omega, arma::mat Sigma1);
arma::vec traceOmegaSigma1(const arma::mat OmegaSigma1, const int numEta);


#endif 
