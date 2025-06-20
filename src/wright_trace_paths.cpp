#include <Rcpp.h>


struct EdgeNum {
  std::string lhs, rhs, op;
  double est;
};


enum class Side : uint8_t { 
  LHS = 0, RHS = 1 
};


inline Side flip(Side s) { 
  return (s == Side::LHS) ? Side::RHS : Side::LHS; 
}


// ─────────────────────────────────────────────────────────────────────────────
//  Tracing coefficients
// ─────────────────────────────────────────────────────────────────────────────
static double trace_rec(const std::string&           x,
                        const std::string&           y,
                        const std::vector<EdgeNum>&  edges,
                        int                          maxLen,
                        double                       curProd,
                        int                          covCount,
                        Side                         side,
                        int                          depth) {
  if (depth > maxLen) {
    Rcpp::warning("Encountered a non-recursive model (infinite loop) when tracing paths");
    return 0.0;
  }

  if (covCount > 1)                return 0.0;       // >=2 covariances -> illegal
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
                     curProd * e.est, // multiply along the path
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

  std::vector<EdgeNum> edges;
  edges.reserve(parTable.nrows());
  for (R_xlen_t i = 0; i < parTable.nrows(); ++i) {
    edges.push_back(EdgeNum { 
      static_cast<std::string>(lhs[i]),
      static_cast<std::string>(rhs[i]),
      static_cast<std::string>(op[i]),
      est[i] 
    });
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


// ─────────────────────────────────────────────────────────────────────────────
//  Tracing labels
// ─────────────────────────────────────────────────────────────────────────────
struct EdgeStr {
  std::string lhs, rhs, op, param;
};


static void trace_rec_paths(const std::string&                        x,
                            const std::string&                        y,
                            const std::vector<EdgeStr>&               edges,
                            int                                       maxLen,
                            std::vector<std::string>&                 curPath,
                            std::vector< std::vector<std::string> >&  outPaths,
                            int                                       covCount,
                            Side                                      side,
                            int                                       depth) {

  if (depth > maxLen) {
    Rcpp::warning("Encountered a non-recursive model (infinite loop) when tracing paths");
    return;
  }
  if (covCount > 1) return;
  if (x == y && side == Side::RHS) {
    outPaths.push_back(curPath);
    return;
  }

  for (const auto& e : edges) {
    bool matches = (side == Side::LHS) ? (e.lhs == x) : (e.rhs == x);
    if (!matches) continue;

    Side        nextSide = side;
    int         nextCov  = covCount;
    std::string nextVar  = (side == Side::LHS) ? e.rhs : e.lhs;

    if (e.op == "~~") {
      nextCov  += 1;
      nextSide  = flip(side);
    }

    curPath.push_back(e.param);
    trace_rec_paths(nextVar, y, edges, maxLen,
                    curPath, outPaths,
                    nextCov, nextSide, depth + 1);
    curPath.pop_back();
  }
}


static Rcpp::List vec2CVList(const std::vector< std::vector<std::string> >& paths) {
  Rcpp::List out(paths.size());
  for (size_t i = 0; i < paths.size(); ++i)
    out[i] = Rcpp::CharacterVector(paths[i].begin(), paths[i].end());
  return out;
}


// [[Rcpp::export]]
Rcpp::List tracePathsCharacterCpp(Rcpp::CharacterVector x,
                                  Rcpp::CharacterVector y,
                                  Rcpp::DataFrame       parTable,
                                  std::string           paramCol = "mod",
                                  int                   maxlen   = 100) {

  if (x.size() != y.size())
    Rcpp::stop("`x` and `y` must have the same length.");

  Rcpp::CharacterVector lhs   = parTable["lhs"];
  Rcpp::CharacterVector rhs   = parTable["rhs"];
  Rcpp::CharacterVector op    = parTable["op"];
  Rcpp::CharacterVector param = parTable[paramCol];

  std::vector<EdgeStr> edges;
  edges.reserve(parTable.nrows());
  for (R_xlen_t i = 0; i < parTable.nrows(); ++i) {
    edges.push_back(EdgeStr { 
        static_cast<std::string>(lhs[i]),
        static_cast<std::string>(rhs[i]),
        static_cast<std::string>(op[i]),
        static_cast<std::string>(param[i]) 
    });
  }

  R_xlen_t n = x.size();
  Rcpp::List out(n);
  Rcpp::CharacterVector nms(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    std::string xi = static_cast<std::string>(x[i]);
    std::string yi = static_cast<std::string>(y[i]);

    std::vector<std::string>                 curPath;
    std::vector< std::vector<std::string> >  foundPaths;

    trace_rec_paths(xi, yi, edges, maxlen,
                    curPath, foundPaths,
                    /*covCount*/ 0,
                    /*side*/     Side::LHS,
                    /*depth*/    1);

    out[i] = vec2CVList(foundPaths);
    nms[i] = xi + "~~" + yi;
  }
  out.attr("names") = nms;
  return out;
}
