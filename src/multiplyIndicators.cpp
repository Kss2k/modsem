#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::NumericVector multiplyIndicatorsCpp(Rcpp::DataFrame df) {
  if (df.size() <= 1) {
    return df[0];
  }

  Rcpp::NumericVector x = df[0];
  Rcpp::NumericVector y = df[1];
  Rcpp::NumericVector product = x*y;

  // Delete the two first entries
  df.erase(df.begin(), df.begin() + 2);
  // Add product to the first entry
  df.push_front(product);

  // Recurse
  return multiplyIndicatorsCpp(df);
}
