#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace std;
using namespace Rcpp;

//' Multiply indicators 
//' @param df A data DataFrame
//' @return A NumericVector
//' @export
// [[Rcpp::export]]
NumericVector multiplyIndicatorsCpp(DataFrame df) {
  if (df.size() <= 1) {
    return df[0];
  }
  NumericVector x = df[0];
  NumericVector y = df[1];
  NumericVector product = x*y;

  // Delete the two first entries
  df.erase(df.begin(), df.begin() + 2);
  // Add product to the first entry
  df.push_front(product);
  // Recursing
  return multiplyIndicatorsCpp(df);
}


