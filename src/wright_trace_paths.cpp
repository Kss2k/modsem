#include <Rcpp.h>


struct Edge {
  std::string lhs, rhs, op;
  double      est;
};

enum class Side : uint8_t { LHS = 0, RHS = 1 };

inline Side flip(Side s) { return (s == Side::LHS) ? Side::RHS : Side::LHS; }


static double trace_rec(const std::string&          x,
                        const std::string&          y,
                        const std::vector<Edge>&    edges,
                        int                         maxLen,
                        double                      curProd,
                        int                         covCount,
                        Side                        side,
                        int                         depth) {
  if (depth > maxLen) {
    Rcpp::warning("Encountered a non-recursive model (infinite loop) when tracing paths");
    return 0.0;
  }
  if (covCount > 1)                return 0.0;         // ≥2 covariances ⇒ illegal
  if (x == y && side == Side::RHS) return curProd;   // finished path

  double sum = 0.0;

  for (const auto& e : edges) {
    bool matches = (side == Side::LHS) ? (e.lhs == x) : (e.rhs == x);
    if (!matches) continue;

    Side        nextSide = side;
    int         nextCov  = covCount;
    std::string nextVar  = (side == Side::LHS) ? e.rhs : e.lhs;

    if (e.op == "~~") { // covariance: flip travel direction
      nextCov  += 1;
      nextSide  = flip(side);
    }

    sum += trace_rec(nextVar, y, edges, maxLen,
                     curProd * e.est,         // multiply along the path
                     nextCov, nextSide, depth + 1);
  }
  return sum;
}


// [[Rcpp::export]]
Rcpp::NumericVector tracePathsNumericCpp(Rcpp::CharacterVector x,
                                         Rcpp::CharacterVector y,
                                         Rcpp::DataFrame       parTable,
                                         int                   maxlen = 100) {
  if (x.size() != y.size())
    Rcpp::stop("`x` and `y` must have the same length.");

  Rcpp::CharacterVector lhs = parTable["lhs"];
  Rcpp::CharacterVector rhs = parTable["rhs"];
  Rcpp::CharacterVector op  = parTable["op"];
  Rcpp::NumericVector   est = parTable["est"];

  std::vector<Edge> edges;
  edges.reserve(parTable.nrows());
  for (R_xlen_t i = 0; i < parTable.nrows(); ++i) {
    edges.push_back( Edge{ static_cast<std::string>(lhs[i]),
                           static_cast<std::string>(rhs[i]),
                           static_cast<std::string>(op[i]),
                           est[i] } );
  }

  R_xlen_t n = x.size();
  Rcpp::NumericVector out(n);
  Rcpp::CharacterVector nms(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    std::string xi = static_cast<std::string>(x[i]);
    std::string yi = static_cast<std::string>(y[i]);

    out[i] = trace_rec(xi, yi, edges, maxlen,
                       /*current product*/ 1.0,
                       /*covCount*/        0,
                       /*side*/            Side::LHS,
                       /*depth*/           1);

    nms[i] = xi + "~~" + yi;  // name:  x~~y
  }

  out.attr("names") = nms;
  return out;
}
