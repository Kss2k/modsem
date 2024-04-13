#ifndef LMS_H
#define LMS_H


arma::vec muLmsCpp(Rcpp::List model, arma::vec z);
arma::mat sigmaLmsCpp(Rcpp::List model, arma::vec z);
arma::mat zToMatrix(arma::vec z, int numEtas);


#endif
