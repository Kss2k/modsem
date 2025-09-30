#ifndef PML_H
#define PML_H


arma::vec fcc_vec_arma(const arma::vec& xj, const arma::vec& xk,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk);

arma::vec foc_vec_arma(const arma::vec& xj, const arma::uvec& r,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk,
                       const arma::vec& tau_k);

arma::vec foo_vec_arma(const arma::uvec& r, const arma::uvec& s,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk,
                       const arma::vec& tau_j, const arma::vec& tau_k);

arma::vec probPML(const arma::mat& data,
                  const arma::vec& mu,
                  const arma::mat& Sigma,
                  const arma::uvec& isOrderedEnum,
                  const arma::mat& thresholds);

#endif
