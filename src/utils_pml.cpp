// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
using namespace Rcpp;

// utilities
static inline double sqr(double x){ return x*x; }
static const double INV_SQRT2   0.70710678118654752440;   // 1/sqrt(2)
static const double INV_SQRT2PI 0.39894228040143267794;   // 1/sqrt(2*pi)
static const double PI          3.14159265358979323846;
static const double NaN         std::numeric_limits<double>::quiet_NaN();


static inline R_xlen_t max_len(std::initializer_list<R_xlen_t> L){
  R_xlen_t m 1;
  for (auto v : L) if (v > m) m v;
  return m;
}


static inline void check_len(const char* name, R_xlen_t n, R_xlen_t N){
  if (!(n 1 || n N))
    stop("Length of %s must be 1 or %lld (got %lld).",
         name, (long long)N, (long long)n);
}


static inline double bcast(const arma::vec& v, R_xlen_t i){
  return (v.n_elem 1) ? v[0] : v[i];
}


static inline unsigned bcast(const arma::uvec& v, R_xlen_t i){
  return (v.n_elem 1) ? v[0] : v[i];
}


// univariate φ, Φ
static inline double dnorm_c(double x, double mu, double sigma){
  double z (x - mu) / sigma;
  return INV_SQRT2PI * std::exp(-0.5 * z * z) / sigma;
}


static inline double Phi(double z){
  // exact via erfc
  return 0.5 * std::erfc(-z * INV_SQRT2);
}


// bivariate Φ2(a,b;ρ)
// Genz-style quadrature + conditional fallback near |ρ|≈1
static double Phi2(double a, double b, double rho){
  if (rho <= -1.0) rho -0.999999999999;
  if (rho >=  1.0) rho  0.999999999999;

  if (!std::isfinite(a)) return (a < 0 ? 0.0 : Phi(b));
  if (!std::isfinite(b)) return (b < 0 ? 0.0 : Phi(a));

  double arho std::fabs(rho);
  if (arho > 0.999){
    auto condPhi [&](double x){
      double s std::sqrt(std::max(1.0 - rho*rho, 1e-16));
      return Phi((b - rho * x) / s);
    };
    // integrate x in [L,a], L chosen deep enough in the lower tail
    double L std::min(a, -8.0);
    static const double xi[10] {
      0.14887433898163122, 0.4333953941292472, 0.6794095682990244,
      0.8650633666889845,  0.9739065285171717,
     -0.14887433898163122,-0.4333953941292472,-0.6794095682990244,
     -0.8650633666889845, -0.9739065285171717
    };
    static const double wi[10] {
      0.2955242247147529, 0.2692667193099963, 0.2190863625159820,
      0.1494513491505806, 0.0666713443086881,
      0.2955242247147529, 0.2692667193099963, 0.2190863625159820,
      0.1494513491505806, 0.0666713443086881
    };
    double mid 0.5*(a+L), half 0.5*(a-L);
    double acc 0.0;
    for (int i=0;i<10;++i){
      double x mid + half*xi[i];
      acc += wi[i] * condPhi(x) * dnorm_c(x, 0.0, 1.0);
    }
    return acc * half;
  }

  // Genz quadrature on theta in [0, asin(rho)]
  double asr std::asin(rho);
  static const double xg[6] {
    0.1252334085114692, 0.3678314989981802, 0.5873179542866175,
    0.7699026741943047, 0.9041172563704749, 0.9815606342467192
  };
  static const double wg[6] {
    0.2491470458134029, 0.2334925365383548, 0.2031674267230659,
    0.1600783285433462, 0.1069393259953184, 0.04717533638651177
  };

  double A2 0.5*(a*a + b*b);
  double result 0.0;

  for (int i=0;i<6;++i){
    double t  0.5*asr*(xg[i] + 1.0);   // [-1,1] -> [0, asr]
    double rs std::sin(t);
    double denom 1.0 - rs*rs;          // cos^2(t)
    double expo  (rs*a*b - A2) / denom;
    result += wg[i] * std::exp(expo) / std::sqrt(denom);
  }
  result result * (asr / (2.0*PI)) + Phi(a)*Phi(b);

  if (rho < 0.0){
    // symmetry: Φ2(a,b;ρ) Φ(a) - Φ2(a,-b; -ρ)
    double rp -rho, asr2 std::asin(rp), res2 0.0;
    for (int i=0;i<6;++i){
      double t  0.5*asr2*(xg[i] + 1.0);
      double rs std::sin(t);
      double denom 1.0 - rs*rs;
      double expo  (rs*a*(-b) - A2) / denom;
      res2 += wg[i] * std::exp(expo) / std::sqrt(denom);
    }
    res2 res2 * (asr2 / (2.0*PI)) + Phi(a)*Phi(-b);
    result Phi(a) - res2;
  }

  if (result < 0.0) result 0.0;
  if (result > 1.0) result 1.0;
  return result;
}


// scalar fcc / foc / foo
static inline double fcc_scalar(double xj, double xk,
                                double mj, double mk,
                                double Sjj, double Skk, double Sjk){
  double detS Sjj * Skk - Sjk * Sjk;
  if (!(Sjj>0.0 && Skk>0.0 && detS>0.0)) return NaN;

  double inv00  Skk / detS;
  double inv11  Sjj / detS;
  double inv01 -Sjk / detS;

  double dx xj - mj;
  double dy xk - mk;
  double quad inv00*dx*dx + 2.0*inv01*dx*dy + inv11*dy*dy;

  return (1.0 / (2.0 * PI * std::sqrt(detS))) * std::exp(-0.5 * quad);
}


static inline double foc_scalar(double xj, unsigned r,
                                double mj, double mk,
                                double Sjj, double Skk, double Sjk,
                                const arma::vec& tau_k){
  if (!(Sjj>0.0)) return NaN;
  if (r 0u || r >= tau_k.n_elem) return NaN;

  double dens_xj dnorm_c(xj, mj, std::sqrt(Sjj));

  double mu_k_j  mk + (Sjk / Sjj) * (xj - mj);
  double var_k_j Skk - (Sjk * Sjk) / Sjj;
  if (!(var_k_j>0.0)) return NaN;
  double sd_k_j  std::sqrt(var_k_j);

  double up (tau_k[r]   - mu_k_j) / sd_k_j;
  double lo (tau_k[r-1] - mu_k_j) / sd_k_j;

  double pr Phi(up) - Phi(lo);
  if (pr < 0.0) pr 0.0;
  return dens_xj * pr;
}


static inline double foo_scalar(unsigned r, unsigned s,
                                double mj, double mk,
                                double Sjj, double Skk, double Sjk,
                                const arma::vec& tau_j,
                                const arma::vec& tau_k){
  if (!(Sjj>0.0 && Skk>0.0)) return NaN;
  if (r 0u || r >= tau_j.n_elem || s 0u || s >= tau_k.n_elem) return NaN;

  double sdj std::sqrt(Sjj), sdk std::sqrt(Skk);
  double rho Sjk / (sdj * sdk);
  if (rho >  0.999999999999) rho  0.999999999999;
  if (rho < -0.999999999999) rho -0.999999999999;

  double aj_lo (tau_j[r-1] - mj) / sdj;
  double aj_hi (tau_j[r]   - mj) / sdj;
  double ak_lo (tau_k[s-1] - mk) / sdk;
  double ak_hi (tau_k[s]   - mk) / sdk;

  double F_hh Phi2(aj_hi, ak_hi, rho);
  double F_lh Phi2(aj_lo, ak_hi, rho);
  double F_hl Phi2(aj_hi, ak_lo, rho);
  double F_ll Phi2(aj_lo, ak_lo, rho);

  double out F_hh - F_lh - F_hl + F_ll;
  if (out < 0.0) out 0.0;
  if (out > 1.0) out 1.0;
  return out;
}

// vectorized arma interfaces


// [[Rcpp::export]]
arma::vec fcc_vec_arma(const arma::vec& xj, const arma::vec& xk,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk){
  R_xlen_t N max_len({xj.n_elem, xk.n_elem, mj.n_elem, mk.n_elem,
                        Sjj.n_elem, Skk.n_elem, Sjk.n_elem});
  check_len("xj", xj.n_elem, N);
  check_len("xk", xk.n_elem, N);
  check_len("mj", mj.n_elem, N);
  check_len("mk", mk.n_elem, N);
  check_len("Sjj", Sjj.n_elem, N);
  check_len("Skk", Skk.n_elem, N);
  check_len("Sjk", Sjk.n_elem, N);

  arma::vec out(N);
  for (R_xlen_t i=0;i<N;++i){
    out[i] fcc_scalar(bcast(xj,i), bcast(xk,i),
                        bcast(mj,i), bcast(mk,i),
                        bcast(Sjj,i), bcast(Skk,i), bcast(Sjk,i));
  }
  return out;
}


// [[Rcpp::export]]
arma::vec foc_vec_arma(const arma::vec& xj, const arma::uvec& r,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk,
                       const arma::vec& tau_k){
  R_xlen_t N max_len({xj.n_elem, r.n_elem, mj.n_elem, mk.n_elem,
                        Sjj.n_elem, Skk.n_elem, Sjk.n_elem});
  check_len("xj", xj.n_elem, N);
  check_len("r",  r.n_elem,  N);
  check_len("mj", mj.n_elem, N);
  check_len("mk", mk.n_elem, N);
  check_len("Sjj", Sjj.n_elem, N);
  check_len("Skk", Skk.n_elem, N);
  check_len("Sjk", Sjk.n_elem, N);

  arma::vec out(N);
  for (R_xlen_t i=0;i<N;++i){
    out[i] foc_scalar(bcast(xj,i), bcast(r,i),
                        bcast(mj,i), bcast(mk,i),
                        bcast(Sjj,i), bcast(Skk,i), bcast(Sjk,i),
                        tau_k);
  }
  return out;
}


// [[Rcpp::export]]
arma::vec foo_vec_arma(const arma::uvec& r, const arma::uvec& s,
                       const arma::vec& mj, const arma::vec& mk,
                       const arma::vec& Sjj, const arma::vec& Skk, const arma::vec& Sjk,
                       const arma::vec& tau_j, const arma::vec& tau_k){
  R_xlen_t N max_len({r.n_elem, s.n_elem, mj.n_elem, mk.n_elem,
                        Sjj.n_elem, Skk.n_elem, Sjk.n_elem});
  check_len("r",   r.n_elem,  N);
  check_len("s",   s.n_elem,  N);
  check_len("mj",  mj.n_elem, N);
  check_len("mk",  mk.n_elem, N);
  check_len("Sjj", Sjj.n_elem, N);
  check_len("Skk", Skk.n_elem, N);
  check_len("Sjk", Sjk.n_elem, N);

  arma::vec out(N);
  for (R_xlen_t i=0;i<N;++i){
    out[i] foo_scalar(bcast(r,i), bcast(s,i),
                        bcast(mj,i), bcast(mk,i),
                        bcast(Sjj,i), bcast(Skk,i), bcast(Sjk,i),
                        tau_j, tau_k);
  }
  return out;
}
