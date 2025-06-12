#ifndef mvnorm_h
#define mvnorm_h 


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);

arma::vec dmvnrm_arma_mc(arma::mat const &x,  arma::rowvec const &mean,  
    arma::mat const &sigma, bool const logd);

double totalDmvnWeightedCpp(const arma::vec& mu, const arma::mat& sigma,
    const arma::vec& nu, const arma::mat& S, double tgamma, int n, int d);

#endif // !mvnorm_h
