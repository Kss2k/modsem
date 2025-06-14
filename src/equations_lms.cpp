#include <RcppArmadillo.h>
#include "lms.h"
#include "mvnorm.h"
#include <float.h>

// [[Rcpp::depends(RcppArmadillo)]]


// Deprecated will remove soon...
// [[Rcpp::export]]
arma::vec muLmsCpp(Rcpp::List model, arma::vec z) {
  Rcpp::List matrices = model["matrices"];
  Rcpp::List info = model["info"];
  Rcpp::List quad = model["quad"];
  int numXis = Rcpp::as<int>(info["numXis"]);
  int k = Rcpp::as<int>(quad["k"]);
  arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  arma::mat Ie = Rcpp::as<arma::mat>(matrices["Ieta"]);
  arma::mat lY = Rcpp::as<arma::mat>(matrices["lambdaY"]);
  arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  arma::mat tY = Rcpp::as<arma::mat>(matrices["tauY"]);
  arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else zVec = arma::zeros<arma::vec>(numXis);
  arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);
  }
  arma::vec muX = tX + lX * (beta0 + A * zVec);
  arma::vec muY = tY + 
    lY * (Binv * (a + 
          Gx * (beta0 + A * zVec) + 
          kronZ.t() * Oxx * (beta0 + A * zVec)));
  return arma::join_cols(muX, muY);
}


// Deprecated will remove soon...
// [[Rcpp::export]]
arma::mat sigmaLmsCpp(Rcpp::List model, arma::vec z) {
  Rcpp::List matrices = model["matrices"];
  Rcpp::List info = model["info"];
  Rcpp::List quad = model["quad"];
  int numXis = Rcpp::as<int>(info["numXis"]);
  int k = Rcpp::as<int>(quad["k"]);
  arma::mat A = Rcpp::as<arma::mat>(matrices["A"]);
  arma::mat Oxx = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
  arma::mat Oex = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
  arma::mat Ie = Rcpp::as<arma::mat>(matrices["Ieta"]);
  arma::mat lY = Rcpp::as<arma::mat>(matrices["lambdaY"]);
  arma::mat lX = Rcpp::as<arma::mat>(matrices["lambdaX"]);
  arma::mat tY = Rcpp::as<arma::mat>(matrices["tauY"]);
  arma::mat tX = Rcpp::as<arma::mat>(matrices["tauX"]);
  arma::mat Gx = Rcpp::as<arma::mat>(matrices["gammaXi"]);
  arma::mat Ge = Rcpp::as<arma::mat>(matrices["gammaEta"]);
  arma::mat a = Rcpp::as<arma::mat>(matrices["alpha"]);
  arma::mat beta0 = Rcpp::as<arma::mat>(matrices["beta0"]);
  arma::mat Psi = Rcpp::as<arma::mat>(matrices["psi"]); 
  arma::mat d = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
  arma::mat e = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);

  arma::vec zVec;
  if (k > 0) zVec = arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else zVec = arma::zeros<arma::vec>(numXis);
  arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

  arma::mat Binv;
  if (Ie.n_cols == 1) {
    Binv = arma::mat(Ie);
  } else {
    Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);
  }

  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));
  arma::mat Sxx = lX * A * Oi * A.t() * lX.t() + d;
  arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
  arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
  arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() + 
    lY * (Binv * Psi * Binv.t()) * lY.t() + e;
  return arma::join_cols(arma::join_rows(Sxx, Sxy), arma::join_rows(Sxy.t(), Syy));
}


inline arma::mat make_Oi(unsigned k, unsigned numXis) {
  arma::mat Oi = arma::eye<arma::mat>(numXis, numXis);
  Oi.diag() = arma::join_cols(arma::zeros<arma::vec>(k), arma::ones<arma::vec>(numXis - k));
  return Oi;
}


inline arma::vec make_zvec(unsigned k, unsigned numXis, const arma::vec& z) {
  if (k > 0) return arma::join_cols(z, arma::zeros<arma::vec>(numXis - k));
  else return arma::zeros<arma::vec>(numXis);
}


struct LMSModel {
  arma::mat A, Oxx, Oex, Ie, lY, lX, tY, tX, Gx, Ge,
    a, beta0, Psi, d, e;
  unsigned  k       = 0;
  unsigned  numXis  = 0;

  explicit LMSModel(const Rcpp::List& modFilled) {

    Rcpp::List matrices = modFilled["matrices"];
    Rcpp::List info     = modFilled["info"];
    Rcpp::List quad     = modFilled["quad"];

    k       = Rcpp::as<unsigned>(quad["k"]);
    numXis  = Rcpp::as<unsigned>(info["numXis"]);

    // one-liners, no loops
    A      = Rcpp::as<arma::mat>(matrices["A"]);
    Oxx    = Rcpp::as<arma::mat>(matrices["omegaXiXi"]);
    Oex    = Rcpp::as<arma::mat>(matrices["omegaEtaXi"]);
    Ie     = Rcpp::as<arma::mat>(matrices["Ieta"]);
    lY     = Rcpp::as<arma::mat>(matrices["lambdaY"]);
    lX     = Rcpp::as<arma::mat>(matrices["lambdaX"]);
    tY     = Rcpp::as<arma::mat>(matrices["tauY"]);
    tX     = Rcpp::as<arma::mat>(matrices["tauX"]);
    Gx     = Rcpp::as<arma::mat>(matrices["gammaXi"]);
    Ge     = Rcpp::as<arma::mat>(matrices["gammaEta"]);
    a      = Rcpp::as<arma::mat>(matrices["alpha"]);
    beta0  = Rcpp::as<arma::mat>(matrices["beta0"]);
    Psi    = Rcpp::as<arma::mat>(matrices["psi"]);
    d      = Rcpp::as<arma::mat>(matrices["thetaDelta"]);
    e      = Rcpp::as<arma::mat>(matrices["thetaEpsilon"]);
  }

  arma::vec mu(const arma::vec& z) const {
    arma::vec zVec = make_zvec(k, numXis, z);
    arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

    arma::mat Binv;
    if (Ie.n_cols == 1) Binv = arma::mat(Ie);
    else Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    arma::vec muX = tX + lX * (beta0 + A * zVec);
    arma::vec muY = tY + 
      lY * (Binv * (a + 
            Gx * (beta0 + A * zVec) + 
            kronZ.t() * Oxx * (beta0 + A * zVec)));
    return arma::join_cols(muX, muY);
  }


  arma::mat Sigma(const arma::vec& z) const {
    arma::vec zVec = make_zvec(k, numXis, z);

    arma::mat kronZ = arma::kron(Ie, beta0 + A * zVec);

    arma::mat Binv;
    if (Ie.n_cols == 1) Binv = arma::mat(Ie);
    else Binv = arma::inv(Ie - Ge - kronZ.t() * Oex);

    arma::mat Oi = make_Oi(k, numXis);
    arma::mat Sxx = lX * A * Oi * A.t() * lX.t() + d;
    arma::mat Eta = Binv * (Gx * A + kronZ.t() * Oxx * A);
    arma::mat Sxy = lX * (A * Oi * Eta.t()) * lY.t();
    arma::mat Syy = lY * Eta * Oi * Eta.t() * lY.t() + 
      lY * (Binv * Psi * Binv.t()) * lY.t() + e;

    return arma::join_cols(
        arma::join_rows(Sxx, Sxy), 
        arma::join_rows(Sxy.t(), Syy)
        );
  }
};


// [[Rcpp::export]]
double completeLogLikLmsCpp(Rcpp::List modelR, Rcpp::List P, Rcpp::List quad) {
  const LMSModel model = LMSModel(modelR);
  const arma::mat Z = Rcpp::as<arma::mat>(quad["n"]).t(); // transpose so we can use column order vectors

  Rcpp::List info   = modelR["info"];
  Rcpp::List TGamma = P["tgamma"];
  Rcpp::List Means  = P["mean"];
  Rcpp::List Covs   = P["cov"];

  const int n = Rcpp::as<int>(info["N"]);
  const int d = Rcpp::as<int>(info["ncol"]);

  double ll = 0;
  for (std::size_t i = 0; i < Z.n_cols; i++) {
    const double tgamma = Rcpp::as<double>(TGamma[i]);
    if (tgamma <= DBL_MIN) continue;

    const arma::vec z = Z.col(i);
    const arma::vec nu = Rcpp::as<arma::vec>(Means[i]);
    const arma::mat S = Rcpp::as<arma::mat>(Covs[i]);
    const arma::mat mu = model.mu(z);
    const arma::mat Sigma = model.Sigma(z);

    ll += totalDmvnWeightedCpp(mu, Sigma, nu, S, tgamma, n, d);
  }

  return ll;
}


struct Jacob { 
  arma::vec dmu; 
  arma::mat dSi; 
};


Jacob jac_beta0(const LMSModel& M,
                const arma::vec& z,          // length k
                std::size_t idx)             // which β₀ entry
{
    Jacob out;

    // --------- dβ --------------------------------------------------
    arma::vec d_z = arma::zeros<arma::vec>(M.beta0.n_elem);
    d_z(idx) = 1.0;

    // --------- μ part ---------------------------------------------
    arma::vec muX = M.lX * d_z;

    arma::vec big_z = make_zvec(M.k, M.numXis, z);
    arma::vec beta  = M.beta0 + M.A * big_z;
    arma::mat kronZ = arma::kron(M.Ie, beta.t());

    arma::mat Binv =
        (M.Ie.n_cols == 1) ? M.Ie
                           : arma::inv(M.Ie - M.Ge - kronZ.t() * M.Oex);

    arma::vec muY = M.lY *
        ( Binv * ( M.Gx * d_z + kronZ.t() * M.Oxx * d_z ) );

    out.dmu = arma::join_cols(muX, muY);

    // --------- Σ part (β₀ does not enter Σ) -----------------------
    out.dSi.zeros(out.dmu.n_elem, out.dmu.n_elem);
    return out;
}


Jacob jac_A(const LMSModel& M,
            const arma::vec& z,        // length k
            std::size_t      r,
            std::size_t      c)
{
    Jacob out;

    arma::vec zBig  = make_zvec(M.k, M.numXis, z);       // ξ vector
    double    zc    = zBig(c);                           // z_c

    // ---------------- dβ (= e_r * z_c) ----------------------------
    arma::vec d_beta(M.beta0.n_elem, arma::fill::zeros);
    d_beta(r) = zc;

    arma::vec  beta   = M.beta0 + M.A * zBig;
    arma::mat  kronZ  = arma::kron(M.Ie, beta);          // (Iη ⊗ β)
    arma::mat  dKronZ = arma::kron(M.Ie, d_beta);        // derivative

    // ---------------- B⁻¹ and its derivative ----------------------
    arma::mat Binv =
        (M.Ie.n_cols == 1)
        ? M.Ie
        : arma::inv(M.Ie - M.Ge - kronZ.t() * M.Oex);

    arma::mat dBinv;
    if (M.Ie.n_cols == 1) {
        dBinv.zeros(Binv.n_rows, Binv.n_cols);           // scalar η
    } else {
        arma::mat dB =              //  dB  =  – (dKronZᵀ Ωηξ)
            - dKronZ.t() * M.Oex;
        dBinv = - Binv * dB * Binv;  //  dB⁻¹ = –B⁻¹ dB B⁻¹   (one minus!)
    }

    // ---------------- μ -part ------------------------------------
    arma::vec muX = M.lX * d_beta;

    arma::vec inner =
        M.a + M.Gx * beta + kronZ.t() * M.Oxx * beta;

    arma::vec dInner =
        M.Gx * d_beta +
        dKronZ.t() * M.Oxx * beta +
        kronZ.t()   * M.Oxx * d_beta;

    arma::vec muY =
        M.lY * ( dBinv * inner + Binv * dInner );

    out.dmu = arma::join_cols(muX, muY);

    // ---------------- Σ -part ------------------------------------
    unsigned p = M.tX.n_elem;
    unsigned q = M.tY.n_elem;
    out.dSi.zeros(p + q, p + q);

    arma::mat Oi   = make_Oi(M.k, M.numXis);

    // ----- dSxx ---------------------------------------------------
    arma::mat Erc(M.A.n_rows, M.A.n_cols, arma::fill::zeros);
    Erc(r, c) = 1.0;

    arma::mat dSxx =
        M.lX * ( Erc * Oi * M.A.t()
               + M.A  * Oi * Erc.t() ) * M.lX.t();

    // ----- Eta and dEta ------------------------------------------
    arma::mat Eta  = Binv *
        ( M.Gx * M.A + kronZ.t() * M.Oxx * M.A );

    arma::mat dEta =
        dBinv * ( M.Gx * M.A + kronZ.t() * M.Oxx * M.A ) +
        Binv  * ( M.Gx * Erc +
                  dKronZ.t() * M.Oxx * M.A +
                  kronZ.t()   * M.Oxx * Erc );

    // ----- dSxy ---------------------------------------------------
    arma::mat dSxy =
        M.lX * ( Erc * Oi * Eta.t()
               + M.A  * Oi * dEta.t() ) * M.lY.t();

    // ----- dSyy ---------------------------------------------------
    arma::mat dSyy =
        M.lY * ( dEta * Oi * Eta.t()
               + Eta  * Oi * dEta.t() ) * M.lY.t() +
        M.lY * ( dBinv * M.Psi * Binv.t()
               + Binv  * M.Psi * dBinv.t() ) * M.lY.t();

    // ----- assemble ------------------------------------------------
    out.dSi.submat(0, 0,     p-1,   p-1)   = dSxx;
    out.dSi.submat(0, p,     p-1,   p+q-1) = dSxy;
    out.dSi.submat(p, 0,     p+q-1, p-1)   = dSxy.t();
    out.dSi.submat(p, p,     p+q-1, p+q-1) = dSyy;

    return out;
}


Jacob jac_lambdaX(const LMSModel& M,
                  const arma::vec& z,
                  std::size_t r,
                  std::size_t c)
{
    Jacob out;
    unsigned p = M.tX.n_elem;
    unsigned q = M.tY.n_elem;

    arma::vec big_z = make_zvec(M.k, M.numXis, z);
    arma::vec beta  = M.beta0 + M.A * big_z;

    // -------- μ part ----------------------------------------------
    out.dmu.zeros(p + q);
    out.dmu(r) = beta(c);                       // only X-part changes

    // -------- Σ part ----------------------------------------------
    arma::mat Oi   = make_Oi(M.k, M.numXis);
    arma::mat Eta  =
        ((M.Ie.n_cols == 1)
         ? M.Ie
         : arma::inv(M.Ie - arma::kron(M.Ie, beta.t()).t() * M.Oex)) *
        ( M.Gx * M.A +
          arma::kron(M.Ie, beta.t()).t() * M.Oxx * M.A );

    arma::mat Erc(p, M.lX.n_cols, arma::fill::zeros);
    Erc(r, c) = 1.0;

    arma::mat dSxx =
        Erc * M.A * Oi * M.A.t() * M.lX.t() +
        M.lX * M.A * Oi * M.A.t() * Erc.t();

    arma::mat dSxy =
        Erc * ( M.A * Oi * Eta.t() ) * M.lY.t();

    out.dSi.zeros(p + q, p + q);
    out.dSi.submat(0, 0,   p-1, p-1)     = dSxx;
    out.dSi.submat(0, p,   p-1, p+q-1)   = dSxy;
    out.dSi.submat(p, 0,   p+q-1, p-1)   = dSxy.t();
    // dSyy = 0 for λX
    return out;
}


Jacob jac_tauX(const LMSModel& M,
               std::size_t r)               // τX_r
{
    Jacob out;
    unsigned p = M.tX.n_elem;
    unsigned q = M.tY.n_elem;

    out.dmu.zeros(p + q);
    out.dmu(r) = 1.0;                        // only that X manifest
    out.dSi.zeros(p + q, p + q);             // Σ unaffected
    return out;
}


Jacob jac_tauY(const LMSModel& M,
               std::size_t r)               // τY_r  (r over Y block!)
{
    Jacob out;
    unsigned p = M.tX.n_elem;
    unsigned q = M.tY.n_elem;

    out.dmu.zeros(p + q);
    out.dmu(p + r) = 1.0;                    // offset by X block
    out.dSi.zeros(p + q, p + q);
    return out;
}


// [[Rcpp::export]]
arma::vec gradLogLikLmsCpp(const Rcpp::List& modelR,
                           const Rcpp::List& P,
                           const arma::uvec& block,
                           const arma::uvec& row,
                           const arma::uvec& col) {
   const LMSModel M = LMSModel(modelR);


    // ---------- unpack EM sufficient statistics -------------------
    Rcpp::List PMean = P["mean"];
    Rcpp::List PCov = P["cov"];
    Rcpp::List PTGamma = P["tgamma"];

    const arma::mat V = Rcpp::as<arma::mat>(P["V"]);
    const std::size_t Mnodes = V.n_rows;
    const std::size_t d = Rcpp::as<arma::vec>(PMean[0]).size();   // manifest dim

    arma::vec g_mu(d, arma::fill::zeros);
    arma::mat G_Si(d, d, arma::fill::zeros);

    for (std::size_t j=0; j<Mnodes; ++j) {

        const arma::vec  nu = Rcpp::as<arma::vec>(PMean[j]);
        const arma::mat  S  = Rcpp::as<arma::mat>(PCov[j]);
        const double     tg = Rcpp::as<double>(PTGamma[j]);

        arma::vec  mu     = M.mu(V.row(j).t());
        arma::mat  Sig    = M.Sigma(V.row(j).t());
        arma::mat  SigInv = arma::inv_sympd(Sig);

        arma::vec r = nu - mu;

        g_mu += tg * SigInv * r;

        arma::mat Delta = S + tg * (r*r.t()) - tg * Sig;

        G_Si += 0.5 * (SigInv * Delta * SigInv);
    }

    // ---------- chain-rule over free parameters -------------------
    std::size_t t = block.n_elem;
    arma::vec grad(t, arma::fill::zeros);

    for (std::size_t k = 0; k < t; ++k) {

        Jacob J;
        switch (block[k]) {
          case 0:  J = jac_beta0   (M, V.row(0).t(), row[k]);                break;
          case 1:  J = jac_A       (M, V.row(0).t(), row[k], col[k]);        break;
          case 2:  J = jac_lambdaX (M, V.row(0).t(), row[k], col[k]);        break;
          case 3:  J = jac_tauX    (M, row[k]);                              break;
          case 4:  J = jac_tauY    (M, row[k]);                              break;
          default: Rcpp::stop("unknown block code");
        }

        grad[k] = arma::dot(g_mu, J.dmu) +
                  arma::accu(G_Si % J.dSi);
    }

    return grad;
}
