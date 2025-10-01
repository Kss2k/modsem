// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <R_ext/RS.h>

using namespace Rcpp;

extern "C" {
  void F77_NAME(pbivnorm)(double* prob,
                          double* lower,
                          double* uppera,
                          double* upperb,
                          int*    infin,
                          double* correl,
                          int*    length);
}

// --------------------- small helpers ---------------------
static inline double Phi_scalar(double z){
  return 0.5 * std::erfc(-z / std::sqrt(2.0));
}
static inline void clamp_rho_vec(arma::vec& r){
  const double eps = 1.0 - 1e-12;
  r.transform([&](double x){
    if (!std::isfinite(x)) return 0.0;
    if (x >  eps) return  eps;
    if (x < -eps) return -eps;
    return x;
  });
}
static inline void check_len(const char* name, R_xlen_t n, R_xlen_t N){
  if (!(n == 1 || n == N))
    stop("Length of %s must be 1 or %lld (got %lld).",
         name, (long long)N, (long long)n);
}
static inline R_xlen_t max_len(std::initializer_list<R_xlen_t> L){
  R_xlen_t m = 1; for (auto v : L) if (v > m) m = v; return m;
}
// recycle/broadcast to length N
static arma::vec expand_to_N(const arma::vec& x, R_xlen_t N){
  if (x.n_elem == N) return x;
  if (x.n_elem == 1) return arma::vec(N).fill(x[0]);
  stop("Incompatible lengths.");
}
static arma::uvec expand_to_N(const arma::uvec& x, R_xlen_t N){
  if (x.n_elem == N) return x;
  if (x.n_elem == 1) return arma::uvec(N).fill(x[0]);
  stop("Incompatible lengths.");
}

// ===================== vectorized arma wrappers (optimized) =====================

// CC: bivariate **density** (no PBIVNORM here). Fully vectorized Armadillo.
// [[Rcpp::export]]
arma::vec fcc_vec_arma(const arma::vec& xj, const arma::vec& xk,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk){
  R_xlen_t N = max_len({xj.n_elem, xk.n_elem, mj.n_elem, mk.n_elem,
                        Sjj.n_elem, Skk.n_elem, Sjk.n_elem});

  check_len("xj",  xj.n_elem,  N);
  check_len("xk",  xk.n_elem,  N);
  check_len("mj",  mj.n_elem,  N);
  check_len("mk",  mk.n_elem,  N);
  check_len("Sjj", Sjj.n_elem, N);
  check_len("Skk", Skk.n_elem, N);
  check_len("Sjk", Sjk.n_elem, N);

  arma::vec XJ  = expand_to_N(xj,  N);
  arma::vec XK  = expand_to_N(xk,  N);
  arma::vec MJ  = expand_to_N(mj,  N);
  arma::vec MK  = expand_to_N(mk,  N);
  arma::vec SJJ = expand_to_N(Sjj, N);
  arma::vec SKK = expand_to_N(Skk, N);
  arma::vec SJK = expand_to_N(Sjk, N);

  arma::vec detS = SJJ % SKK - SJK % SJK;
  arma::vec inv00 = SKK / detS;
  arma::vec inv11 = SJJ / detS;
  arma::vec inv01 = -SJK / detS;

  arma::vec dx = XJ - MJ;
  arma::vec dy = XK - MK;
  arma::vec quad = inv00 % (dx % dx) + 2.0 * inv01 % (dx % dy) + inv11 % (dy % dy);

  arma::vec out(N, arma::fill::zeros);
  arma::uvec ok = arma::find( (SJJ > 0) % (SKK > 0) % (detS > 0) );
  if (!ok.is_empty()) {
    arma::vec det_ok = detS.elem(ok);
    arma::vec q_ok   = quad.elem(ok);
    out.elem(ok) = (1.0 / (2.0 * M_PI)) * arma::exp(-0.5 * q_ok) / arma::sqrt(det_ok);
  }
  out.elem(arma::find_nonfinite(out)).zeros(); // NaN -> 0
  return out;
}

// OC/CO: density(continuous) * Pr(interval | continuous) ; vectorized.
// [[Rcpp::export]]
arma::vec foc_vec_arma(const arma::vec& xj, const arma::uvec& r,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk,
                       const arma::vec& tau_k){
  R_xlen_t N = max_len({xj.n_elem, r.n_elem, mj.n_elem, mk.n_elem,
                        Sjj.n_elem, Skk.n_elem, Sjk.n_elem});

  check_len("xj",  xj.n_elem,  N);
  check_len("r",   r.n_elem,   N);
  check_len("mj",  mj.n_elem,  N);
  check_len("mk",  mk.n_elem,  N);
  check_len("Sjj", Sjj.n_elem, N);
  check_len("Skk", Skk.n_elem, N);
  check_len("Sjk", Sjk.n_elem, N);

  arma::vec XJ  = expand_to_N(xj,  N);
  arma::uvec R  = expand_to_N(r,   N);
  arma::vec MJ  = expand_to_N(mj,  N);
  arma::vec MK  = expand_to_N(mk,  N);
  arma::vec SJJ = expand_to_N(Sjj, N);
  arma::vec SKK = expand_to_N(Skk, N);
  arma::vec SJK = expand_to_N(Sjk, N);

  arma::vec sdj = arma::sqrt(SJJ);
  arma::vec dens_xj = (1.0 / (std::sqrt(2.0 * M_PI))) * arma::exp(-0.5 * arma::square((XJ - MJ) / sdj)) / sdj;

  arma::vec mu_k_j  = MK + (SJK / SJJ) % (XJ - MJ);
  arma::vec var_k_j = SKK - (SJK % SJK) / SJJ;

  arma::vec out(N, arma::fill::zeros);
  arma::uvec ok = arma::find(var_k_j > 0.0);
  if (!ok.is_empty()){
    arma::vec sd_k_j = arma::sqrt(var_k_j.elem(ok));
    arma::vec up(ok.n_elem), lo(ok.n_elem);

    // Map category index r_t (1..K) to thresholds lo=tau[r-1], up=tau[r]
    for (arma::uword i = 0; i < ok.n_elem; ++i){
      arma::uword t = ok[i];
      unsigned rr = R[t];
      if (rr == 0u || rr >= tau_k.n_elem) { up[i] = NAN; lo[i] = NAN; continue; }
      up[i] = (tau_k[rr]   - mu_k_j[t]) / sd_k_j[i];
      lo[i] = (tau_k[rr-1] - mu_k_j[t]) / sd_k_j[i];
    }

    // Phi(up) - Phi(lo)
    arma::vec pr(ok.n_elem, arma::fill::zeros);
    for (arma::uword i = 0; i < ok.n_elem; ++i){
      if (std::isfinite(up[i]) && std::isfinite(lo[i])) {
        pr[i] = Phi_scalar(up[i]) - Phi_scalar(lo[i]);
        if (pr[i] < 0.0) pr[i] = 0.0;
      }
    }
    out.elem(ok) = dens_xj.elem(ok) % pr;
  }
  return out;
}

// OO: rectangle probability via ONE batched PBIVNORM call (length = 4N).
// [[Rcpp::export]]
arma::vec foo_vec_arma(const arma::uvec& r, const arma::uvec& s,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk,
                       const arma::vec& tau_j, const arma::vec& tau_k){
  R_xlen_t N = max_len({r.n_elem, s.n_elem, mj.n_elem, mk.n_elem,
                        Sjj.n_elem, Skk.n_elem, Sjk.n_elem});
  check_len("r",   r.n_elem,  N);
  check_len("s",   s.n_elem,  N);
  check_len("mj",  mj.n_elem, N);
  check_len("mk",  mk.n_elem, N);
  check_len("Sjj", Sjj.n_elem, N);
  check_len("Skk", Skk.n_elem, N);
  check_len("Sjk", Sjk.n_elem, N);

  arma::uvec R  = expand_to_N(r,   N);
  arma::uvec S  = expand_to_N(s,   N);
  arma::vec MJ  = expand_to_N(mj,  N);
  arma::vec MK  = expand_to_N(mk,  N);
  arma::vec SJJ = expand_to_N(Sjj, N);
  arma::vec SKK = expand_to_N(Skk, N);
  arma::vec SJK = expand_to_N(Sjk, N);

  arma::vec sdj = arma::sqrt(SJJ);
  arma::vec sdk = arma::sqrt(SKK);
  arma::vec rho = SJK / (sdj % sdk);
  clamp_rho_vec(rho);

  // Build 4N corner arrays (order: hh, lh, hl, ll for each t)
  const int len = static_cast<int>(4 * N);
  std::vector<double> UPA(len), UPB(len), COR(len), PROB(len, 0.0);

  for (arma::uword t = 0; t < N; ++t){
    unsigned rr = R[t], ss = S[t];
    // invalid categories -> probability 0
    if (rr == 0u || rr >= tau_j.n_elem || ss == 0u || ss >= tau_k.n_elem) {
      UPA[4*t+0] = UPA[4*t+1] = UPA[4*t+2] = UPA[4*t+3] = 0.0;
      UPB[4*t+0] = UPB[4*t+1] = UPB[4*t+2] = UPB[4*t+3] = 0.0;
      COR[4*t+0] = COR[4*t+1] = COR[4*t+2] = COR[4*t+3] = 0.0;
      continue;
    }

    double aj_lo = (tau_j[rr-1] - MJ[t]) / sdj[t];
    double aj_hi = (tau_j[rr]   - MJ[t]) / sdj[t];
    double ak_lo = (tau_k[ss-1] - MK[t]) / sdk[t];
    double ak_hi = (tau_k[ss]   - MK[t]) / sdk[t];
    double rtt   = rho[t];

    // (hh, lh, hl, ll)
    UPA[4*t+0] = aj_hi; UPB[4*t+0] = ak_hi; COR[4*t+0] = rtt;
    UPA[4*t+1] = aj_lo; UPB[4*t+1] = ak_hi; COR[4*t+1] = rtt;
    UPA[4*t+2] = aj_hi; UPB[4*t+2] = ak_lo; COR[4*t+2] = rtt;
    UPA[4*t+3] = aj_lo; UPB[4*t+3] = ak_lo; COR[4*t+3] = rtt;
  }

  // Call PBIVNORM once
  double LOWER[2] = {0.0, 0.0};  // ignored for INFIN=0
  int    INFIN[2] = {0, 0};
  int    LEN      = len;
  F77_NAME(pbivnorm)(PROB.data(), LOWER, UPA.data(), UPB.data(), INFIN, COR.data(), &LEN);

  // Combine by inclusion–exclusion
  arma::vec out(N, arma::fill::zeros);
  for (arma::uword t = 0; t < N; ++t){
    double val = PROB[4*t+0] - PROB[4*t+1] - PROB[4*t+2] + PROB[4*t+3];
    if (std::isfinite(val)) {
      if (val < 0.0) val = 0.0;
      if (val > 1.0) val = 1.0;
      out[t] = val;
    } else {
      out[t] = 0.0;
    }
  }
  return out;
}


// ---- vectorized safe log (replaces per-element loop) -----------------
static inline arma::vec safe_log(const arma::vec& v, double eps = 1e-300){
  arma::vec out(v.n_elem);
  out.fill(std::log(eps));
  arma::uvec ok = arma::find(v > 0.0);     // (NaN, <=0) -> log(eps)
  if (!ok.is_empty()) out.elem(ok) = arma::log(v.elem(ok));
  return out;
}


// [[Rcpp::export]]
arma::vec probPML(
    const arma::mat& data,                 // (N x P)
    const arma::vec& mu,                   // length P
    const arma::mat& Sigma,                // (P x P)
    const arma::uvec& isOrderedEnum,       // length P; 0=cont, >0=1-based thr idx
    const arma::mat& thresholds
){
  const arma::uword N = data.n_rows;
  const arma::uword P = data.n_cols;

  arma::vec out_log(N, arma::fill::zeros);

  // Split columns once
  std::vector<arma::vec>  ccols(P); // continuous cols (as double)
  std::vector<arma::uvec> ocols(P); // ordinal cols (as 1..K indices)
  for (arma::uword j = 0; j < P; ++j) {
    if (isOrderedEnum[j]) {
      ocols[j] = arma::conv_to<arma::uvec>::from(data.col(j));
    } else {
      ccols[j] = data.col(j);
    }
  }

  // Pairwise accumulation (vectorized by row)
  for (arma::uword j = 1; j < P; ++j) {
    const bool oj = (isOrderedEnum[j] != 0u);
    const double mj = mu[j];
    const double Sjj = Sigma(j,j);

    for (arma::uword i = 0; i < j; ++i) {
      const bool oi = (isOrderedEnum[i] != 0u);

      // Pair constants as length-1 vecs (broadcasted inside *_vec_arma)
      arma::vec mi(1);  mi[0]  = mu[i];
      arma::vec mjv(1); mjv[0] = mj;
      arma::vec Sii(1); Sii[0] = Sigma(i,i);
      arma::vec Sjjv(1);Sjjv[0]= Sjj;
      arma::vec Sij(1); Sij[0] = Sigma(i,j);

      arma::vec term_prob; term_prob.zeros(N);

      if (oi && oj) {
        // OO: one batched PBIVNORM call per pair via foo_vec_arma
        const arma::uvec& r = ocols[i];
        const arma::uvec& s = ocols[j];
        const arma::vec& tau_i = thresholds.row(isOrderedEnum[i]-1).t();
        const arma::vec& tau_j = thresholds.row(isOrderedEnum[j]-1).t();

        // Rcpp::Rcout << "tau_i:\n" << tau_i << "\n";
        // Rcpp::Rcout << "tau_j:\n" << tau_j << "\n";
        // Rcpp::Rcout << "Sii:\n" << Sii << "\n";
        // Rcpp::Rcout << "Sij:\n" << Sij << "\n";
        // Rcpp::Rcout << "Sjjv:\n" << Sjjv << "\n";
        // Rcpp::Rcout << "r:\n" << r << "\n";
        // Rcpp::Rcout << "s:\n" << s << "\n";
        term_prob = foo_vec_arma(r, s, mi, mjv, Sii, Sjjv, Sij, tau_i, tau_j);

      } else if (oi && !oj) {
        // OC: density(x_j) * Pr(interval | x_j)
        const arma::uvec& r = ocols[i];
        const arma::vec&   x = ccols[j];
        const arma::vec& tau_i = thresholds.row(isOrderedEnum[i]-1).t();
        // foc expects (x_cont, r_ord, mj, mk, Sjj, Skk, Sjk, tau_k)
        term_prob = foc_vec_arma(x, r, mjv, mi, Sjjv, Sii, Sij, tau_i);

      } else if (!oi && oj) {
        // CO: symmetric to above
        const arma::vec&   x = ccols[i];
        const arma::uvec& s  = ocols[j];
        const arma::vec& tau_j = thresholds.row(isOrderedEnum[j]-1).t();
        term_prob = foc_vec_arma(x, s, mi, mjv, Sii, Sjjv, Sij, tau_j);

      } else {
        // CC: bivariate normal **density** (already vectorized)
        const arma::vec& x = ccols[i];
        const arma::vec& y = ccols[j];
        term_prob = fcc_vec_arma(x, y, mi, mjv, Sii, Sjjv, Sij);
      }

      // Vectorized accumulation in log-domain
      out_log += safe_log(term_prob);
    }
  }

  return out_log;
}



