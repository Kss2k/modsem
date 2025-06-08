// -----------------------------------------------------------------------------
//  lms_aghq.cpp – Rcpp/Armadillo back-end for adaptive-quadrature LMS
//  (June 2025, patched for numerical robustness)
// -----------------------------------------------------------------------------
#include <RcppArmadillo.h>
#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::mat;  using arma::vec;  using arma::uword;  using arma::span;

/*───────────────────────────────────────────────────────────────────────────*
 | 0.  Robust helpers                                                       |
 *───────────────────────────────────────────────────────────────────────────*/
inline mat safe_inv_sympd(const mat& M, double ridge = 1e-8)
{
  mat out;
  if (!arma::inv_sympd(out, M)) {                       // 1st try
    mat Mp = M + ridge * arma::eye(M.n_rows, M.n_cols);
    if (!arma::inv_sympd(out, Mp))                      // 2nd try
      out = arma::pinv(Mp);                             // pseudo-inverse
  }
  return out;
}

inline mat safe_chol(const mat& M,
                     double ridge = 1e-8, unsigned int maxit = 6)
{
  mat Mp = M;
  for (unsigned int t = 0; t < maxit; ++t) {
    try { return arma::chol(Mp); }                      // success → return
    catch (...) {
      double bump = ridge * std::pow(10.0, static_cast<int>(t));
      Mp = M + bump * arma::eye(M.n_rows, M.n_cols);
    }
  }
  /* eigen fallback */
  vec  d;  mat Q;
  arma::eig_sym(d, Q, M);
  d.transform([&](double v){ return v < ridge ? ridge : v; });
  mat pd = Q * arma::diagmat(d) * Q.t();
  return arma::chol(pd);
}

/*───────────────────────────────────────────────────────────────────────────*
 | 1.  Fast MVN density                                                     |
 *───────────────────────────────────────────────────────────────────────────*/
inline double dmvnorm_fast(const vec& x, const vec& mu, const mat& Sigma)
{
  static const double log2pi = std::log(2.0 * M_PI);
  uword d = x.n_elem;
  double sign_det, log_det;
  arma::log_det(log_det, sign_det, Sigma);
  vec diff = x - mu;
  double quad = arma::as_scalar(diff.t() * safe_inv_sympd(Sigma) * diff);
  return std::exp(-0.5 * (d*log2pi + log_det + quad));
}

/*───────────────────────────────────────────────────────────────────────────*
 | 2.  Immutable filled-model container                                     |
 *───────────────────────────────────────────────────────────────────────────*/
class LmsModel {
public:
  /* structural matrices */
  mat A, omegaXiXi, omegaEtaXi, Ieta,
      lambdaY, lambdaX, tauY, tauX,
      gammaXi, gammaEta, alpha, psi,
      thetaDelta, thetaEpsilon;
  uword k, numXis;

  LmsModel(const Rcpp::List& mod)
  {
    Rcpp::List M = mod["matrices"];
    #define PULL(X)  X = Rcpp::as<mat>(M[#X]);
    PULL(A)  PULL(omegaXiXi) PULL(omegaEtaXi) PULL(Ieta)
    PULL(lambdaY) PULL(lambdaX) PULL(tauY) PULL(tauX)
    PULL(gammaXi) PULL(gammaEta) PULL(alpha) PULL(psi)
    PULL(thetaDelta) PULL(thetaEpsilon)
    #undef PULL

    k = Rcpp::as<Rcpp::List>(mod["quad"]).containsElementNamed("k")
          ? static_cast<uword>(Rcpp::as<int>(Rcpp::List(mod["quad"])["k"])) : 1;
    numXis = Rcpp::as<Rcpp::List>(mod["info"]).containsElementNamed("numXis")
          ? static_cast<uword>(Rcpp::as<int>(Rcpp::List(mod["info"])["numXis"])) : k;
  }

  /* μ(z₁) */
  vec mu(const vec& z1) const
  {
    vec zVec(numXis, arma::fill::zeros);                // pad zeros
    zVec(span(0, k-1)) = z1;

    mat kronZ = kron(Ieta, A * zVec);
    mat Binv  = (Ieta.n_cols == 1) ? Ieta
      : safe_inv_sympd(Ieta - gammaEta - kronZ.t() * omegaEtaXi);

    vec muX = tauX + lambdaX * A * zVec;
    vec muY = tauY +
      lambdaY * (Binv * (alpha + gammaXi * A * zVec + kronZ.t()*omegaXiXi*A*zVec));

    return arma::join_vert(muX, muY);
  }

  /* Σ(z₁) */
  mat sigma(const vec& z1) const
  {
    vec zVec(numXis, arma::fill::zeros);
    zVec(span(0, k-1)) = z1;

    mat kronZ = kron(Ieta, A * zVec);
    mat Binv  = (Ieta.n_cols == 1) ? Ieta
      : safe_inv_sympd(Ieta - gammaEta - kronZ.t() * omegaEtaXi);

    mat OI = arma::eye(numXis, numXis);
    for (uword i = 0; i < k; ++i) OI(i,i) = 0.0;

    mat Sxx = lambdaX * A * OI * A.t() * lambdaX.t() + thetaDelta;
    mat G   = Binv * (gammaXi * A + kronZ.t() * omegaXiXi * A);
    mat Sxy = lambdaX * A * OI * G.t() * lambdaY.t();
    mat Syy = lambdaY * (G*OI*G.t() + Binv*psi*Binv.t()) * lambdaY.t() + thetaEpsilon;

    mat out(Sxx.n_rows + Syy.n_rows, Sxx.n_rows + Syy.n_rows, arma::fill::zeros);
    out(span(0,Sxx.n_rows-1), span(0,Sxx.n_cols-1)) = Sxx;
    out(span(0,Sxy.n_rows-1), span(Sxx.n_cols,out.n_cols-1)) = Sxy;
    out(span(Sxx.n_rows,out.n_rows-1), span(0,Sxy.n_rows-1)) = Sxy.t();
    out(span(Sxx.n_rows,out.n_rows-1), span(Sxx.n_cols,out.n_cols-1)) = Syy;
    return out;
  }
};

/*───────────────────────────────────────────────────────────────────────────*
 | 3.  Gradient / Hessian for mode-finding                                   |
 *───────────────────────────────────────────────────────────────────────────*/
struct LogF {
  const LmsModel& mod; const arma::rowvec& y;
  LogF(const LmsModel& m, const arma::rowvec& yr): mod(m), y(yr) {}
  double operator()(const vec& z) const {
    return std::log(dmvnorm_fast(vec(y.t()), mod.mu(z), mod.sigma(z))) +
           arma::accu(arma::log_normpdf(z));
  }
};

void grad_hess(const LogF& f, const vec& z, vec& g, mat& H, double eps=1e-5)
{
  uword k = z.n_elem;
  g.set_size(k);  H.set_size(k,k);
  double f0 = f(z);
  for (uword i=0;i<k;++i){
    vec z_u=z; z_u(i)+=eps;  vec z_d=z; z_d(i)-=eps;
    double fu=f(z_u), fd=f(z_d);
    g(i)=(fu-fd)/(2*eps);
    for(uword j=i;j<k;++j){
      vec z_pp=z_u; z_pp(j)+=eps;
      vec z_mm=z_d; z_mm(j)-=eps;
      double fpp=f(z_pp);
      H(i,j) = (fpp-fu-fd+f0)/(eps*eps);
      H(j,i) = H(i,j);
    }
  }
}

/*───────────────────────────────────────────────────────────────────────────*
 | 4.  Standard tensor-product GH grid                                       |
 *───────────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
Rcpp::List gh_rule_cpp(unsigned int k, unsigned int m=5)
{
  vec nodes(m), w1(m);
  for (unsigned int j=0;j<m;++j){
    double x = std::sqrt(2.0)*R::qnorm((j+0.75)/(m+0.5),0,1,1,0);
    nodes(j)=x;  w1(j)=std::sqrt(M_PI)/(m+0.5);
  }
  unsigned int M = std::pow(m,k);
  mat S(M,k); vec w(M,arma::fill::ones);
  for (unsigned int idx=0;idx<M;++idx){
    unsigned int r=idx;
    for(unsigned int d=0;d<k;++d){
      unsigned int pos=r % m;
      S(idx,d)=nodes(pos);  w(idx)*=w1(pos);  r/=m;
    }
  }
  return Rcpp::List::create(Rcpp::Named("S")=S, Rcpp::Named("w")=w);
}

/*───────────────────────────────────────────────────────────────────────────*
 | 5.  Single-case adaptive grid                                             |
 *───────────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
Rcpp::List adapt_case_cpp(const Rcpp::List& modelR, const arma::rowvec& y,
                          const mat& S, const vec& w_std,
                          unsigned int maxit=40)
{
  LmsModel mod(modelR);
  uword k=S.n_cols;
  vec z(k,arma::fill::zeros);
  LogF f(mod,y);

  /* Newton */
  for(unsigned int it=0;it<maxit;++it){
    vec g; mat H; grad_hess(f,z,g,H);
    vec step = arma::solve(-H, g, arma::solve_opts::likely_sympd);
    if(arma::norm(step,2) < 1e-6) break;
    z += step;
  }
  vec g; mat H; grad_hess(f,z,g,H);
  mat L = safe_chol(safe_inv_sympd(H)/2.0);

  mat Z = S*L.t(); Z.each_row() += z.t();
  vec w = w_std * arma::det(L);
  return Rcpp::List::create(Rcpp::Named("Z")=Z,Rcpp::Named("w")=w);
}

/*───────────────────────────────────────────────────────────────────────────*
 | 6.  QuadCache & XPtr                                                      |
 *───────────────────────────────────────────────────────────────────────────*/
struct QuadCache {
  std::vector<mat> Zs; std::vector<vec> ws, Ps; std::vector<double> ll;
  QuadCache(size_t n): Zs(n),ws(n),Ps(n),ll(n) {}
};
using quadPtr = Rcpp::XPtr<QuadCache>;

/*───────────────────────────────────────────────────────────────────────────*
 | 7.  E-step: build cache                                                   |
 *───────────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
SEXP estep_cpp(const Rcpp::List& modelR, const arma::mat& data,
               const mat& S, const vec& w_std)
{
  LmsModel mod(modelR);  size_t N=data.n_rows;
  quadPtr cache(new QuadCache(N), true);

  for(size_t i=0;i<N;++i){
    auto ad = adapt_case_cpp(modelR, data.row(i), S, w_std);
    mat Z = ad["Z"]; vec w = ad["w"];
    size_t m_i=w.n_elem; vec p(m_i);
    for(size_t j=0;j<m_i;++j)
      p(j)=dmvnorm_fast(vec(data.row(i).t()), mod.mu(Z.row(j).t()),
                        mod.sigma(Z.row(j).t())) * w(j);
    double fi = arma::accu(p);
    cache->ll[i]=fi;  p/=fi;
    cache->Zs[i]=Z;  cache->ws[i]=w;  cache->Ps[i]=p;
  }
  return cache;
}

/*───────────────────────────────────────────────────────────────────────────*
 | 8.  Complete-data log-lik for M-step                                      |
 *───────────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
double loglik_mstep_cpp(const Rcpp::List& modelR, SEXP ptr)
{
  quadPtr c(ptr);  LmsModel mod(modelR); size_t N=c->Zs.size(); double ll=0;
  for(size_t i=0;i<N;++i){
    const mat& Z=c->Zs[i]; const vec& p=c->Ps[i];
    vec contrib(p.n_elem);
    for(size_t j=0;j<p.n_elem;++j){
      contrib(j)=std::log(dmvnorm_fast(vec(Z.row(j).t()),
                                       mod.mu(Z.row(j).t()),
                                       mod.sigma(Z.row(j).t())))+std::log(p(j));
    }
    double m=contrib.max(); ll += m+std::log(arma::accu(arma::exp(contrib-m)));
  }
  return ll;
}

/*───────────────────────────────────────────────────────────────────────────*
 | 9.  Cache utilities                                                       |
 *───────────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]] double obs_loglik_from_cache(SEXP ptr){
  quadPtr c(ptr); return std::accumulate(c->ll.begin(), c->ll.end(), 0.0,
      [](double a,double fi){return a+std::log(fi);} );}

// [[Rcpp::export]] void free_cache(SEXP ptr){ quadPtr(ptr).release(); }

/*───────────────────────────────────────────────────────────────────────────*
 |10.  Observed LL reusing grid but updated θ                                |
 *───────────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
double obs_loglik_from_grid_cpp(const Rcpp::List& modelR, SEXP ptr)
{
  quadPtr c(ptr);  LmsModel mod(modelR); size_t N=c->Zs.size(); double ll=0;
  for(size_t i=0;i<N;++i){
    const mat& Z=c->Zs[i]; const vec& w=c->ws[i];
    vec dens(w.n_elem);
    for(size_t j=0;j<w.n_elem;++j)
      dens(j)=dmvnorm_fast(vec(Z.row(j).t()), mod.mu(Z.row(j).t()),
                           mod.sigma(Z.row(j).t()));
    ll += std::log(arma::dot(w,dens));
  }
  return ll;
}

/*───────────────────────────────────────────────────────────────────────────*
 |11.  Differentiable observed LL (rebuilds grid each call)                  |
 *───────────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
double obs_loglik_cpp(const Rcpp::List& modelR, const arma::mat& data,
                      const mat& S, const vec& w_std, unsigned int maxit=40)
{
  LmsModel mod(modelR); size_t N=data.n_rows; double ll=0;
  for(size_t i=0;i<N;++i){
    auto ad=adapt_case_cpp(modelR,data.row(i),S,w_std,maxit);
    mat Z=ad[\"Z\"]; vec w=ad[\"w\"]; vec dens(w.n_elem);
    for(size_t j=0;j<w.n_elem;++j)
      dens(j)=dmvnorm_fast(vec(data.row(i).t()), mod.mu(Z.row(j).t()),
                           mod.sigma(Z.row(j).t()));
    ll+=std::log(arma::dot(w,dens));
  }
  return ll;
}
// ---------------------------------------------------------------------------
//  (end of file)
// ---------------------------------------------------------------------------