#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "qml.h"
#include "mvnorm.h"

// [[Rcpp::export]]
arma::mat muQmlCpp(Rcpp::List m, int t, int ncores = 1) {
  ThreadSetter ts(ncores);                       // set threads

  const int  numEta      = m["numEta"];
  const int  numXi       = m["numXi"];
  const arma::mat alpha  = m["alpha"];
  const arma::mat beta0  = m["beta0"];
  const arma::mat gammaXi= m["gammaXi"];
  const arma::mat omegaXiXi = m["omegaXiXi"];
  const arma::mat L1     = m["L1"];
  const arma::mat L2     = m["L2"];
  const arma::mat X      = m["x"];
  const arma::mat U      = m["u"];
  const arma::mat Sigma1 = m["Sigma1"];
  const arma::mat Binv   = m["Binv"];
  const arma::mat kronXi = m["kronXi"];

  const arma::vec trOmegaSigma = traceOmegaSigma1(omegaXiXi * Sigma1, numEta);
  arma::mat Ey(t, numEta, arma::fill::none);

  const int lastColKOxx = numXi * numEta - 1,
            lastColBinv = numEta - 1;

  if (Binv.n_rows > static_cast<unsigned>(numEta)) {
    #ifdef _OPENMP
    #pragma omp parallel for if(ncores>1) schedule(static)
    #endif

    for (int i = 0; i < t; ++i) {
      const int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;

      const arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOxx);
      const arma::mat Binv_t   = Binv  .submat(firstRow, 0, lastRow, lastColBinv);

      Ey.row(i) = ( Binv_t *
        ( trOmegaSigma + alpha
          + gammaXi * (beta0 + L1 * X.row(i).t())
          + kronXi_t * omegaXiXi * (beta0 + L1 * X.row(i).t()) )
        + L2 * U.row(i).t()
      ).t();
    }

  } else {
    #ifdef _OPENMP
    #pragma omp parallel for if(ncores>1) schedule(static)
    #endif

    for (int i = 0; i < t; ++i) {
      const int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;
      const arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOxx);

      Ey.row(i) = ( Binv *
        ( trOmegaSigma + alpha
          + gammaXi * (beta0 + L1 * X.row(i).t())
          + kronXi_t * omegaXiXi * (beta0 + L1 * X.row(i).t()) )
        + L2 * U.row(i).t()
      ).t();
    }
  }

  return Ey;
}


// [[Rcpp::export]]
arma::mat sigmaQmlCpp(Rcpp::List m, int t, int ncores = 1) {
  ThreadSetter ts(ncores);

  const int  numEta           = m["numEta"];
  const int  numXi            = m["numXi"];
  const arma::mat gammaXi     = m["gammaXi"];
  const arma::mat omegaXiXi   = m["omegaXiXi"];
  const arma::mat L1          = m["L1"];
  const arma::mat L2          = m["L2"];
  const arma::mat X           = m["x"];
  const arma::mat U           = m["u"];
  const arma::mat Sigma1      = m["Sigma1"];
  const arma::mat Sigma2Theta = m["Sigma2ThetaEpsilon"];
  const arma::mat psi         = m["psi"];
  const arma::mat Binv        = m["Binv"];
  const arma::mat kronXi      = m["kronXi"];

  const arma::mat varZ = varZCpp(omegaXiXi, Sigma1, numEta);
  arma::mat sigmaE(t * numEta, numEta, arma::fill::none);

  const int lastColSigma = numEta - 1,
            lastColKOxx  = numXi * numEta - 1;

  const arma::mat omegaXiXi2T = omegaXiXi + transposeOmega(omegaXiXi, numEta); // omega + omega'

  if (Binv.n_rows > static_cast<unsigned>(numEta)) {
    #ifdef _OPENMP
    #pragma omp parallel for if(ncores>1) schedule(static)
    #endif

    for (int i = 0; i < t; ++i) {
      const int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;
      const arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOxx);
      const arma::mat Binv_t   = Binv  .submat(firstRow, 0, lastRow, lastColSigma);
      const arma::mat Sigma2   = Binv_t * psi * Binv_t.t() + Sigma2Theta;

      const arma::mat BinvGammaXi2Omega =
        Binv_t * gammaXi + kronXi_t * omegaXiXi2T;

      sigmaE.submat(firstRow, 0, lastRow, lastColSigma) =
        BinvGammaXi2Omega * Sigma1 * BinvGammaXi2Omega.t() +
        Sigma2 + Binv_t * varZ * Binv_t.t();
    }

  } else {
    const arma::mat BinvVarZ = Binv * varZ * Binv.t();
    const arma::mat Sigma2   = Binv * psi * Binv.t() + Sigma2Theta;

    #ifdef _OPENMP
    #pragma omp parallel for if(ncores>1) schedule(static)
    #endif

    for (int i = 0; i < t; ++i) {
      const int firstRow = i * numEta, lastRow = (i + 1) * numEta - 1;
      const arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOxx);
      const arma::mat BinvGammaXi2Omega =
        Binv * gammaXi + kronXi_t * omegaXiXi2T;

      sigmaE.submat(firstRow, 0, lastRow, lastColSigma) =
        BinvGammaXi2Omega * Sigma1 * BinvGammaXi2Omega.t() +
        Sigma2 + BinvVarZ;
    }
  }

  return sigmaE;
}


// [[Rcpp::export]]
arma::mat calcKronXi(Rcpp::List m, int t, int ncores = 1) {
  ThreadSetter ts(ncores);

  const int numEta  = m["numEta"];
  const int numXi   = m["numXi"];
  const arma::mat beta0 = m["beta0"];
  const arma::mat L1    = m["L1"];
  const arma::mat X     = m["x"];
  const arma::mat Ie    = m["Ieta"];

  arma::mat out(t * numEta, numXi * numEta, arma::fill::none);

  #ifdef _OPENMP
  #pragma omp parallel for if(ncores>1) schedule(static)
  #endif

  for (int i = 0; i < t; ++i) {
    out.submat(i * numEta, 0,
               (i + 1) * numEta - 1, numXi * numEta - 1) =
      arma::kron(Ie, beta0.t() + X.row(i) * L1.t());
  }

  return out;
}


// [[Rcpp::export]]
arma::mat calcBinvCpp(Rcpp::List m, int t, int ncores = 1) {
  ThreadSetter ts(ncores);

  const int numEta = m["numEta"],
            numXi  = m["numXi"],
            kOmega = m["kOmegaEta"];

  const arma::mat gammaEta   = m["gammaEta"];
  const arma::mat Ie         = m["Ieta"];
  const arma::mat B0         = Ie - gammaEta; // Invariant
  const arma::mat omegaEtaXi = m["omegaEtaXi"];

  if (numEta == 1) return Ie;
  if (kOmega == 0) return arma::inv(B0);

  const arma::mat kronXi = m["kronXi"];
  arma::mat Bout(t * numEta, numEta, arma::fill::none);

  const int lastColB   = numEta - 1,
            lastColKOx = numXi * numEta - 1;

  #ifdef _OPENMP
  #pragma omp parallel for if(ncores>1) schedule(static)
  #endif

  for (int i = 0; i < t; ++i) {
    const int firstRow = i * numEta,
              lastRow  = (i + 1) * numEta - 1;

    const arma::mat kronXi_t = kronXi.submat(firstRow, 0, lastRow, lastColKOx);

    Bout.submat(firstRow, 0, lastRow, lastColB) =
      arma::inv(B0 - kronXi_t * omegaEtaXi);
  }

  return Bout;
}


// [[Rcpp::export]]
arma::vec dnormCpp(const arma::vec& x, const arma::vec& mu, const arma::vec& sigma,
                   int ncores = 1) {
  ThreadSetter ts(ncores);

  const double log2pi = std::log(2.0 * M_PI);
  arma::vec out(x.n_rows, arma::fill::none);

  #ifdef _OPENMP
  #pragma omp parallel for if(ncores>1) schedule(static)
  #endif

  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const double diff = x(i) - mu(i);
    const double sd   = sigma(i);

    out(i) = -0.5 * log2pi - std::log(sd) - 0.5 * (diff * diff) / (sd * sd);
  }

  return out;
}


// [[Rcpp::export]]
arma::mat varZCpp(arma::mat Omega, arma::mat Sigma1, int numEta) {
  arma::mat varZ = arma::mat(numEta, numEta);
  const int subRows = Omega.n_rows / numEta;

  for (int i = 0; i < numEta; i++) {
    varZ(i, i) = varZSubOmega(Omega.submat(i * subRows, 0,
          (i + 1) * subRows - 1, (Omega.n_cols - 1)), Sigma1);
  }

  return varZ;
}


double varZSubOmega(arma::mat Omega, arma::mat Sigma1) {
  const int ds = Sigma1.n_rows;
  double varZ = 0;

  for (int i = 0; i < ds; i++) {
    for (int j = 0; j < ds; j++) {
      for (int k = 0; k < ds; k++) {
        for (int s = 0; s < ds; s++) {
          varZ += Omega(i, j) * Omega(k, s) *
            (Sigma1(i, j) * Sigma1(k, s) +
             Sigma1(i, k) * Sigma1(j, s) +
             Sigma1(i, s) * Sigma1(j, k));
        }
      }
    }
  }

  const double trOmegaSigma1 = arma::trace(Omega * Sigma1);
  return varZ - trOmegaSigma1 * trOmegaSigma1;
}


arma::vec traceOmegaSigma1(const arma::mat OmegaSigma1, const int numEta) {
  arma::vec   trace = arma::vec(numEta);
  const int subRows = OmegaSigma1.n_rows / numEta;

  for (int i = 0; i < numEta; i++) {
    for (int j = 0; j < int(OmegaSigma1.n_cols); j++) {
      trace(i) += OmegaSigma1(i * subRows + j, j);
    }
  }

  return trace;
}


arma::mat transposeOmega(const arma::mat Omega, const int numEta) {
  if (numEta <= 1) return Omega.t();

  const int subRows = Omega.n_rows / numEta;
  arma::mat OmegaT(Omega.n_rows, Omega.n_cols, arma::fill::zeros);

  for (int i = 0; i < numEta; i++) {
    const int ri = i * subRows;
    const int ci = 0;
    const int rk = (i + 1) * subRows - 1;
    const int ck = Omega.n_cols - 1;

    OmegaT.submat(ri, ci, rk, ck) = Omega.submat(ri, ci, rk, ck).t();
  }

  return OmegaT;
}


// =============================================================
// QMLModel: parameter struct with precomputed cache
// =============================================================

inline arma::uvec qmlCharMatch(const Rcpp::CharacterVector& needle,
                               const Rcpp::CharacterVector& haystack) {
  arma::uvec idx(needle.size());
  for (int k = 0; k < needle.size(); ++k) {
    bool found = false;
    for (int j = 0; j < haystack.size(); ++j) {
      if (needle[k] == haystack[j]) { idx(k) = (arma::uword)j; found = true; break; }
    }
    if (!found)
      Rcpp::stop("qmlCharMatch: '%s' not found",
                 Rcpp::as<std::string>(needle[k]).c_str());
  }
  return idx;
}

inline Rcpp::CharacterVector qmlMatColnames(SEXP x) {
  Rcpp::RObject ro(x);
  SEXP dn = ro.attr("dimnames");
  return Rcpp::as<Rcpp::CharacterVector>(VECTOR_ELT(dn, 1));
}


struct QMLModel {
  // Free parameter matrices  (DA_BLOCKS: 0=lX,1=lY,2=tX,3=tY,4=d,5=e,
  //   7=Psi,8=a,9=beta0,10=Gx,11=Ge,12=Oxx,13=Oex,16=phi)
  arma::mat lX, lY, tX, tY, d, e, Psi, a, beta0, Gx, Ge, Oxx, Oex, phi, Ie;

  unsigned numXis, numEtas;
  int      kOmegaEta;
  bool     hasR;   // any latent etas with multiple indicators
  bool     hasTE2; // any non-latent etas needing residual variance in fullSig2TE

  // Column indices into data.full
  arma::uvec x_col_idx;
  arma::uvec y_col_idx;
  arma::vec  sampling_weights;
  bool       has_sw;

  // Structural indices (hasR)
  arma::uvec R_na_linidx;    // col-major NA positions in emptyR
  arma::uvec lY_free_linidx; // col-major !selectScalingY positions in lY
  arma::uvec colsR_uvec;     // positions of colsR in allIndsEtas (= e diag idx = y col idx)
  arma::uvec betaRows_uvec;  // row positions of selectBetaRows in lY
  arma::uvec latentEta_uvec; // positions of latentEtas in etas ordering
  arma::uvec selectTE1_diag; // diagonal positions in e for scaling-indicator residuals

  // Structural indices (hasTE2)
  arma::uvec selectTE2_diag;    // diagonal positions in e for non-latent scaling residuals
  arma::uvec nonLatentEta_uvec; // positions of non-latent etas in etas ordering

  // Indices for f3: nonNormalInds = allIndsEtas NOT in colsU (scaling/single indicators)
  arma::uvec nonNormal_uvec;    // positions within etas_cv (= y column ordering)

  // Template for R (0 where fixed, will be overwritten at NA positions)
  arma::mat R_tmpl;

  // Cached obs-independent derived quantities
  arma::vec  tauX_adj;
  arma::mat  LXPLX, cholLXPLX, invLXPLX, L1, Sigma1;    // cholLXPLX is upper-triangular
  arma::mat  R, RER, cholRER, invRER;           // cholRER is upper-triangular
  arma::mat  Beta, subTE1, L2_sub, Sigma2TE_sub;
  arma::mat  fullSig2TE;                        // numEtas x numEtas
  arma::mat  L2R_cache;                         // numEtas x nColsR
  arma::vec  trOmSig;
  arma::mat  varZ, Oxx2T;
  arma::mat  Binv_c, Sigma2_c, BinvVarZ_c;     // kOmegaEta == 0 only

  void update_cache() {
    tauX_adj = tX + lX * beta0;

    LXPLX     = lX * phi * lX.t() + d;
    cholLXPLX = arma::chol(LXPLX, "upper");
    invLXPLX  = arma::inv_sympd(LXPLX);
    L1        = phi * lX.t() * invLXPLX;
    Sigma1    = phi - L1 * lX * phi;

    if (hasR) {
      R = R_tmpl;
      for (arma::uword k = 0; k < R_na_linidx.n_elem; ++k)
        R(R_na_linidx(k)) = -lY(lY_free_linidx(k));

      RER    = R * e.submat(colsR_uvec, colsR_uvec) * R.t();
      cholRER = arma::chol(RER, "upper");
      invRER = arma::inv_sympd(RER);

      Beta = lY.submat(betaRows_uvec, latentEta_uvec);

      const arma::uword ns = selectTE1_diag.n_elem;
      arma::vec te1(ns);
      for (arma::uword k = 0; k < ns; ++k)
        te1(k) = e(selectTE1_diag(k), selectTE1_diag(k));
      subTE1 = arma::diagmat(te1);

      L2_sub       = -subTE1 * Beta.t() * invRER;
      Sigma2TE_sub = subTE1 - arma::diagmat(arma::square(te1)) * Beta.t() * invRER * Beta;

      fullSig2TE.zeros();
      fullSig2TE.submat(latentEta_uvec, latentEta_uvec) = Sigma2TE_sub;
      L2R_cache.zeros();
      L2R_cache.rows(latentEta_uvec) = L2_sub * R;
    }

    if (hasTE2) {
      for (arma::uword k = 0; k < selectTE2_diag.n_elem; ++k) {
        const arma::uword p = nonLatentEta_uvec(k);
        fullSig2TE(p, p) = e(selectTE2_diag(k), selectTE2_diag(k));
      }
    }

    trOmSig = traceOmegaSigma1(Oxx * Sigma1, numEtas);
    varZ    = varZCpp(Oxx, Sigma1, numEtas);
    Oxx2T   = Oxx + transposeOmega(Oxx, numEtas);

    if (kOmegaEta == 0) {
      Binv_c     = arma::inv(Ie - Ge);
      Sigma2_c   = Binv_c * Psi * Binv_c.t() + fullSig2TE;
      BinvVarZ_c = Binv_c * varZ * Binv_c.t();
    }
  }

  explicit QMLModel(const Rcpp::List& submodel) {
    const Rcpp::List mat  = submodel["matrices"];
    const Rcpp::List info = submodel["info"];
    const Rcpp::List dat  = submodel["data"];

    numXis    = Rcpp::as<unsigned>(info["numXis"]);
    numEtas   = Rcpp::as<unsigned>(info["numEtas"]);
    kOmegaEta = Rcpp::as<int>(info["kOmegaEta"]);

    lX    = Rcpp::as<arma::mat>(mat["lambdaX"]);
    lY    = Rcpp::as<arma::mat>(mat["lambdaY"]);
    tX    = Rcpp::as<arma::mat>(mat["tauX"]);
    tY    = Rcpp::as<arma::mat>(mat["tauY"]);
    d     = Rcpp::as<arma::mat>(mat["thetaDelta"]);
    e     = Rcpp::as<arma::mat>(mat["thetaEpsilon"]);
    Psi   = Rcpp::as<arma::mat>(mat["psi"]);
    a     = Rcpp::as<arma::mat>(mat["alpha"]);
    beta0 = Rcpp::as<arma::mat>(mat["beta0"]);
    Gx    = Rcpp::as<arma::mat>(mat["gammaXi"]);
    Ge    = Rcpp::as<arma::mat>(mat["gammaEta"]);
    Oxx   = Rcpp::as<arma::mat>(mat["omegaXiXi"]);
    Oex   = Rcpp::as<arma::mat>(mat["omegaEtaXi"]);
    phi   = Rcpp::as<arma::mat>(mat["phi"]);
    Ie    = Rcpp::as<arma::mat>(mat["Ieta"]);

    const Rcpp::CharacterVector xis_cv  = info["allIndsXis"];
    const Rcpp::CharacterVector etas_cv = info["allIndsEtas"];
    const Rcpp::CharacterVector eta_lv  = info["etas"];  // latent variable names

    const Rcpp::CharacterVector data_cols = qmlMatColnames(dat["data.full"]);
    x_col_idx = qmlCharMatch(xis_cv,  data_cols);
    y_col_idx = qmlCharMatch(etas_cv, data_cols);

    // nonNormal_uvec: positions in etas_cv NOT absorbed into u (i.e. not in colsU).
    // These are the scaling / single-indicator etas used for the f3 density.
    {
      const SEXP colsU_s = mat["colsU"];
      if (!Rf_isNull(colsU_s) && Rf_length(colsU_s) > 0) {
        const Rcpp::CharacterVector colsU_cv =
          Rcpp::as<Rcpp::CharacterVector>(colsU_s);
        arma::uvec is_u(etas_cv.size(), arma::fill::zeros);
        for (int j = 0; j < etas_cv.size(); ++j)
          for (int k = 0; k < colsU_cv.size(); ++k)
            if (etas_cv[j] == colsU_cv[k]) { is_u(j) = 1; break; }
        nonNormal_uvec = arma::find(is_u == 0);
      } else {
        nonNormal_uvec =
          arma::regspace<arma::uvec>(0, (arma::uword)etas_cv.size() - 1);
      }
    }

    has_sw = !Rf_isNull(dat["weights"]);
    if (has_sw) sampling_weights = Rcpp::as<arma::vec>(dat["weights"]);

    fullSig2TE.zeros(numEtas, numEtas);

    const SEXP emptyR_s = mat["emptyR"];
    hasR = !Rf_isNull(emptyR_s);

    if (hasR) {
      const arma::mat emptyR_mat = Rcpp::as<arma::mat>(emptyR_s);
      R_na_linidx = arma::find_nonfinite(emptyR_mat);
      R_tmpl = emptyR_mat;
      R_tmpl.replace(arma::datum::nan, 0.0);

      lY_free_linidx = arma::find(
          Rcpp::as<arma::mat>(mat["selectScalingY"]) < 0.5);

      colsR_uvec =
          qmlCharMatch(Rcpp::as<Rcpp::CharacterVector>(mat["colsR"]), etas_cv);

      betaRows_uvec = arma::find(
          Rcpp::as<arma::mat>(mat["selectBetaRows"]).col(0) > 0.5);

      const Rcpp::CharacterVector latEta_cv = info["latentEtas"];
      latentEta_uvec = qmlCharMatch(latEta_cv, eta_lv);

      selectTE1_diag = arma::find(
          Rcpp::as<arma::mat>(mat["selectThetaEpsilon1"]).diag() > 0.5);

      L2R_cache.zeros(numEtas, colsR_uvec.n_elem);
    }

    hasTE2 = false;
    const SEXP selTE2_s = mat["selectThetaEpsilon2"];
    if (!Rf_isNull(selTE2_s)) {
      const arma::mat selTE2 = Rcpp::as<arma::mat>(selTE2_s);
      if (selTE2.n_elem > 0) {
        selectTE2_diag = arma::find(selTE2.diag() > 0.5);
        hasTE2 = (selectTE2_diag.n_elem > 0);
      }
    }

    if (hasTE2) {
      if (hasR) {
        arma::uvec isLat(numEtas, arma::fill::zeros);
        for (arma::uword k = 0; k < latentEta_uvec.n_elem; ++k)
          isLat(latentEta_uvec(k)) = 1;
        nonLatentEta_uvec = arma::find(isLat == 0);
      } else {
        nonLatentEta_uvec = arma::regspace<arma::uvec>(0, numEtas - 1);
      }
    }

    update_cache();
  }

  QMLModel thread_clone() const {
    QMLModel c = *this;
    c.lX = arma::mat(lX);  c.lY = arma::mat(lY);
    c.tX = arma::mat(tX);  c.tY = arma::mat(tY);
    c.d  = arma::mat(d);   c.e  = arma::mat(e);
    c.Psi   = arma::mat(Psi);   c.a     = arma::mat(a);
    c.beta0 = arma::mat(beta0); c.Gx    = arma::mat(Gx);
    c.Ge    = arma::mat(Ge);    c.Oxx   = arma::mat(Oxx);
    c.Oex   = arma::mat(Oex);   c.phi   = arma::mat(phi);
    return c;
  }
};


inline double& qml_param(QMLModel& M, std::size_t blk,
                         std::size_t r, std::size_t c) {
  switch (blk) {
    case 0 : return M.lX   (r,c);
    case 1 : return M.lY   (r,c);
    case 2 : return M.tX   (r,c);
    case 3 : return M.tY   (r,c);
    case 4 : return M.d    (r,c);
    case 5 : return M.e    (r,c);
    case 7 : return M.Psi  (r,c);
    case 8 : return M.a    (r,c);
    case 9 : return M.beta0(r,c);
    case 10: return M.Gx   (r,c);
    case 11: return M.Ge   (r,c);
    case 12: return M.Oxx  (r,c);
    case 13: return M.Oex  (r,c);
    case 16: return M.phi  (r,c);
    default: Rcpp::stop("qml_param: unknown block %zu", blk);
  }
}


inline double logLikQmlFromModel(const QMLModel& M,
                                 const arma::mat& data_full,
                                 const int ncores = 1) {
  const int t      = static_cast<int>(data_full.n_rows);
  const int pX     = static_cast<int>(M.x_col_idx.n_elem);
  const int pY     = static_cast<int>(M.y_col_idx.n_elem);
  const int pY_f3  = static_cast<int>(M.nonNormal_uvec.n_elem); // scaling/single inds
  static const double log2pi_d = std::log(2.0 * M_PI);

  arma::mat X(t, pX), Y(t, pY);
  for (int j = 0; j < pX; ++j)
    X.col(j) = data_full.col(M.x_col_idx(j)) - M.tauX_adj(j);
  for (int j = 0; j < pY; ++j)
    Y.col(j) = data_full.col(M.y_col_idx(j)) - M.tY(j, 0);

  // Pre-compute constant parts of f2 log-densities
  const double logdet_x  = 2.0 * arma::sum(arma::log(M.cholLXPLX.diag()));
  const double f2x_const = -0.5 * (pX * log2pi_d + logdet_x);

  int    nU = 0;
  double f2u_const = 0.0;
  arma::mat Y_colsR;
  if (M.hasR) {
    nU = static_cast<int>(M.RER.n_rows);
    const double logdet_u = 2.0 * arma::sum(arma::log(M.cholRER.diag()));
    f2u_const = -0.5 * (nU * log2pi_d + logdet_u);

    const arma::uword nColsR = M.colsR_uvec.n_elem;
    Y_colsR.set_size(t, nColsR);
    for (arma::uword j = 0; j < nColsR; ++j)
      Y_colsR.col(j) = Y.col(M.colsR_uvec(j));
  }

  const bool kOmega0 = (M.kOmegaEta == 0);
  double ll   = 0.0;
  bool   fail = false;

#ifdef _OPENMP
  const bool in_parallel = omp_in_parallel();
  const int  nth = in_parallel ? 1 : ncores;
#else
  const int  nth = 1;
#endif

#pragma omp parallel for if(nth > 1) schedule(static) \
  reduction(+:ll) reduction(||:fail)
  for (int i = 0; i < t; ++i) {
    const arma::vec xi = X.row(i).t();
    // yi: only scaling/single indicators (nonNormalInds), one per eta
    arma::vec yi(pY_f3);
    for (int jj = 0; jj < pY_f3; ++jj) yi(jj) = Y(i, M.nonNormal_uvec(jj));

    // --- f2x: logN(xi | 0, LXPLX) ---
    // cholLXPLX is upper triangular U: LXPLX = U^T U
    // maha = ||U^{-T} xi||^2, solved via trimatl(U^T)
    const arma::vec zx = arma::solve(arma::trimatl(M.cholLXPLX.t()), xi,
                                     arma::solve_opts::fast);
    const double f2x_i = f2x_const - 0.5 * arma::dot(zx, zx);

    // --- f2u: logN(u_raw | 0, RER) ---
    double f2u_i = 0.0;
    if (M.hasR) {
      const arma::vec u_raw = M.R * Y_colsR.row(i).t();
      const arma::vec zu    = arma::solve(arma::trimatl(M.cholRER.t()), u_raw,
                                          arma::solve_opts::fast);
      f2u_i = f2u_const - 0.5 * arma::dot(zu, zu);
    }

    // --- f3: logN(yi | Ey_i, SigmaE_i) ---
    const arma::vec mi = M.beta0 + M.L1 * xi;
    const arma::mat Ki = arma::kron(M.Ie, mi.t());

    arma::mat Binv_i, Sig2_i, BvZ_i;
    if (kOmega0) {
      Binv_i = M.Binv_c;
      Sig2_i = M.Sigma2_c;
      BvZ_i  = M.BinvVarZ_c;
    } else {
      Binv_i = arma::inv(M.Ie - M.Ge - Ki * M.Oex);
      Sig2_i = Binv_i * M.Psi * Binv_i.t() + M.fullSig2TE;
      BvZ_i  = Binv_i * M.varZ * Binv_i.t();
    }

    arma::vec Ey_i = Binv_i * (M.trOmSig + M.a + M.Gx * mi + Ki * M.Oxx * mi);
    if (M.hasR)
      Ey_i += M.L2R_cache * Y_colsR.row(i).t();

    const arma::mat BG2O = Binv_i * M.Gx + Ki * M.Oxx2T;
    const arma::mat SE_i = BG2O * M.Sigma1 * BG2O.t() + Sig2_i + BvZ_i;

    arma::mat Lf;
    if (!arma::chol(Lf, SE_i, "lower")) { fail = true; continue; }

    const arma::vec dv   = yi - Ey_i;
    const arma::vec zf   = arma::solve(arma::trimatl(Lf), dv, arma::solve_opts::fast);
    const double logdet_f = 2.0 * arma::sum(arma::log(Lf.diag()));
    const double f3_i    = -0.5 * (pY_f3 * log2pi_d + logdet_f + arma::dot(zf, zf));

    const double w = M.has_sw ? M.sampling_weights(i) : 1.0;
    ll += w * (f2x_i + f2u_i + f3_i);
  }

  if (fail) return arma::datum::nan;
  return ll;
}


// [[Rcpp::export]]
double logLikQmlCpp(const Rcpp::List& submodel, const int ncores = 1) {
  try {
    const QMLModel M(submodel);
    const arma::mat data_full =
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);
    return logLikQmlFromModel(M, data_full, ncores);
  } catch (...) {
    // Return NaN so the optimizer can step away from this bad parameter point.
    // Do NOT call Rcpp::stop() — that would kill the optimizer via tryCatch.
    return arma::datum::nan;
  }
}


// =============================================================
// Parameter get/set helpers
// =============================================================

inline arma::vec getParamsQml(QMLModel& M,
                               const arma::uvec& block,
                               const arma::uvec& row,
                               const arma::uvec& col) {
  const std::size_t p = block.n_elem;
  arma::vec pars(p);
  for (std::size_t k = 0; k < p; ++k)
    pars[k] = qml_param(M, block[k], row[k], col[k]);
  return pars;
}

inline void setParamsQml(QMLModel& M,
                          const arma::uvec& block,
                          const arma::uvec& row,
                          const arma::uvec& col,
                          const arma::uvec& symmetric,
                          const arma::vec& vals) {
  const std::size_t p = block.n_elem;
  for (std::size_t k = 0; k < p; ++k) {
    qml_param(M, block[k], row[k], col[k]) = vals[k];
    if (symmetric[k] && row[k] != col[k])
      qml_param(M, block[k], col[k], row[k]) = vals[k];
  }
}


// =============================================================
// FD gradient fully in C++ (parallelize over parameters)
// =============================================================

// Inner core: no ThreadSetter, so safe to call from inside OMP regions.
// try-catch inside the loop prevents Armadillo exceptions from escaping OMP.
template<class F>
arma::vec gradientFDQmlCore(QMLModel& M, F logLik,
    const arma::uvec& block, const arma::uvec& row, const arma::uvec& col,
    const arma::uvec& symmetric, const double eps, const int ncores) {
  const std::size_t p = block.n_elem;
  arma::vec grad(p, arma::fill::zeros);
  double f0 = arma::datum::nan;
  try { f0 = logLik(M); } catch (...) {}

  #pragma omp parallel for if(ncores > 1) schedule(static) \
    shared(M, block, row, col, symmetric, eps, grad, f0, p) firstprivate(logLik)
  for (std::size_t k = 0; k < p; ++k) {
    try {
      QMLModel Mc = M.thread_clone();
      double& ti = qml_param(Mc, block[k], row[k], col[k]);
      double* tj = nullptr;
      if (symmetric[k] && row[k] != col[k])
        tj = &qml_param(Mc, block[k], col[k], row[k]);
      ti += eps;
      if (tj) *tj += eps;
      Mc.update_cache();
      grad[k] = (logLik(Mc) - f0) / eps;
    } catch (...) {
      grad[k] = arma::datum::nan;
    }
  }
  return grad;
}


// =============================================================
// Analytical gradient: single-pass accumulator
// Steps 3-5 will add the parameter mappings below.
// =============================================================

arma::vec analyticalGradQmlCore(
    const QMLModel& M,
    const arma::mat& data_full,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric)
{
  const int        t      = static_cast<int>(data_full.n_rows);
  const int        pX     = static_cast<int>(M.x_col_idx.n_elem);
  const int        pY     = static_cast<int>(M.y_col_idx.n_elem);
  const int        pY_f3  = static_cast<int>(M.nonNormal_uvec.n_elem);
  const std::size_t p     = block.n_elem;
  const bool       kOmega0 = (M.kOmegaEta == 0);

  arma::vec grad(p, arma::fill::zeros);

  // Mean-centered indicators
  arma::mat X(t, pX), Y(t, pY);
  for (int j = 0; j < pX; ++j)
    X.col(j) = data_full.col(M.x_col_idx(j)) - M.tauX_adj(j);
  for (int j = 0; j < pY; ++j)
    Y.col(j) = data_full.col(M.y_col_idx(j)) - M.tY(j, 0);

  const arma::uword nColsR = M.hasR ? M.colsR_uvec.n_elem : 0;
  arma::mat Y_colsR;
  if (M.hasR) {
    Y_colsR.set_size(t, nColsR);
    for (arma::uword j = 0; j < nColsR; ++j)
      Y_colsR.col(j) = Y.col(M.colsR_uvec(j));
  }

  // -----------------------------------------------------------------------
  // Accumulators
  // -----------------------------------------------------------------------
  double    sumW  = 0.0;
  arma::mat sumXX(pX, pX, arma::fill::zeros);
  arma::mat sumUU, sumYY;
  arma::vec sumY_colsR_vec;
  if (M.hasR) {
    sumUU.zeros(M.RER.n_rows, M.RER.n_rows);
    sumYY.zeros(M.colsR_uvec.n_elem, M.colsR_uvec.n_elem);
    sumY_colsR_vec.zeros(nColsR);
  }

  // f3 score accumulators
  arma::vec s_Ey_acc(pY_f3,        arma::fill::zeros);
  arma::mat S_SE_acc(pY_f3, pY_f3, arma::fill::zeros);
  // Step 5b: score w.r.t. Sigma1, accumulated as BG2O_i' · S_SE_i · BG2O_i
  arma::mat S_SE_BG2O_acc(M.Gx.n_cols, M.Gx.n_cols, arma::fill::zeros);

  bool fail = false;

  for (int i = 0; i < t; ++i) {
    const arma::vec xi = X.row(i).t();
    const double    w  = M.has_sw ? M.sampling_weights(i) : 1.0;
    sumW  += w;
    sumXX += w * (xi * xi.t());

    if (M.hasR) {
      const arma::vec yc    = Y_colsR.row(i).t();
      const arma::vec u_raw = M.R * yc;
      sumUU += w * (u_raw * u_raw.t());
      sumYY += w * (yc * yc.t());
      sumY_colsR_vec += w * yc;
    }

    // f3 intermediates — same computation as logLikQmlFromModel
    arma::vec yi(pY_f3);
    for (int jj = 0; jj < pY_f3; ++jj) yi(jj) = Y(i, M.nonNormal_uvec(jj));

    const arma::vec mi = M.beta0 + M.L1 * xi;
    const arma::mat Ki = arma::kron(M.Ie, mi.t());

    arma::mat Binv_i, Sig2_i, BvZ_i;
    if (kOmega0) {
      Binv_i = M.Binv_c; Sig2_i = M.Sigma2_c; BvZ_i = M.BinvVarZ_c;
    } else {
      Binv_i = arma::inv(M.Ie - M.Ge - Ki * M.Oex);
      Sig2_i = Binv_i * M.Psi * Binv_i.t() + M.fullSig2TE;
      BvZ_i  = Binv_i * M.varZ * Binv_i.t();
    }

    arma::vec Ey_i = Binv_i * (M.trOmSig + M.a + M.Gx * mi + Ki * M.Oxx * mi);
    if (M.hasR) Ey_i += M.L2R_cache * Y_colsR.row(i).t();

    const arma::mat BG2O = Binv_i * M.Gx + Ki * M.Oxx2T;
    const arma::mat SE_i = BG2O * M.Sigma1 * BG2O.t() + Sig2_i + BvZ_i;

    arma::mat invSE_i;
    if (!arma::inv_sympd(invSE_i, SE_i)) { fail = true; continue; }

    const arma::vec dv_i    = yi - Ey_i;
    const arma::mat S_SE_i  = (-0.5) * (invSE_i - invSE_i * dv_i * dv_i.t() * invSE_i);
    s_Ey_acc      += w * (invSE_i * dv_i);
    S_SE_acc      += w * S_SE_i;
    S_SE_BG2O_acc += w * (BG2O.t() * S_SE_i * BG2O);
  }

  if (fail) { grad.fill(arma::datum::nan); return grad; }

  // -----------------------------------------------------------------------
  // f2x: score w.r.t. LXPLX
  // -----------------------------------------------------------------------
  const arma::mat S_LXPLX = -0.5 * sumW * M.invLXPLX
                           + 0.5 * M.invLXPLX * sumXX * M.invLXPLX;

  // Step 3: map S_LXPLX → grad entries for lX(0), d(4), phi(16)
  // lX(r,c):   dLL/d(lX_{r,c})   = 2 * (S_LXPLX * lX * phi)(r,c)
  // d(r,r):    dLL/d(d_{r,r})    = S_LXPLX(r,r)
  // phi(r,c):  dLL/d(phi_{r,c})  = fac * (lX' * S_LXPLX * lX)(r,c),
  //            fac = 2 when off-diagonal symmetric, 1 otherwise
  {
    const arma::mat SlX_phi  = S_LXPLX * M.lX * M.phi;    // pX  × pXi
    const arma::mat lXT_S_lX = M.lX.t() * S_LXPLX * M.lX; // pXi × pXi

    for (std::size_t k = 0; k < p; ++k) {
      const arma::uword r   = row[k];
      const arma::uword c   = col[k];
      const bool        sym = static_cast<bool>(symmetric[k]);
      const double      fac = (sym && r != c) ? 2.0 : 1.0;

      switch (block[k]) {
        case 0:  grad[k] += 2.0 * SlX_phi(r, c);    break; // lX
        case 4:  grad[k] += S_LXPLX(r, c);           break; // d
        case 16: grad[k] += fac * lXT_S_lX(r, c);   break; // phi
        default: break;
      }
    }
  }

  // -----------------------------------------------------------------------
  // f2u: score w.r.t. RER
  // -----------------------------------------------------------------------
  if (M.hasR) {
    const arma::mat S_RER = -0.5 * sumW * M.invRER
                          + 0.5 * M.invRER * sumUU * M.invRER;

    // Step 4: map S_RER → grad entries for lY(1) and e(5) in colsR.
    //
    // RER = R * e_sub * R',  e_sub = e[colsR_uvec, colsR_uvec].
    //
    // e(q,q) at colsR position rr:
    //   d(RER)/d(e_{q,q}) = R * E_{rr,rr} * R'
    //   grad += (R' * S_RER * R)(rr, rr)
    //
    // lY(r,c) free — maps to R(Ri,Rj) = -lY(r,c):
    //   d(RER)/d(lY_{r,c}) = -(E_{Ri,Rj}*e_sub*R' + R*e_sub*E_{Rj,Ri})
    //   grad += -2 * (S_RER * R * e_sub)(Ri, Rj)
    //   If symmetric: also add the (c,r) partner's contribution.
    {
      const arma::mat e_sub        = M.e.submat(M.colsR_uvec, M.colsR_uvec);
      const arma::mat SRR_e        = S_RER * M.R * e_sub;         // nU × nColsR
      const arma::mat invRER_R_YY  = M.invRER * M.R * sumYY;      // nU × nColsR
      const arma::mat RT_S_R       = M.R.t() * S_RER * M.R;       // nColsR × nColsR

      // Lookup: lY linear col-major index → R linear col-major index.
      // Unused entries stay at sentinel (arma::uword)-1.
      const arma::uword lY_sz = M.lY.n_rows * M.lY.n_cols;
      std::vector<arma::uword> lY2R(lY_sz, (arma::uword)-1);
      for (arma::uword kk = 0; kk < M.lY_free_linidx.n_elem; ++kk)
        lY2R[M.lY_free_linidx(kk)] = M.R_na_linidx(kk);

      const arma::uword nColsR = M.colsR_uvec.n_elem;
      const arma::uword nU     = M.R.n_rows;

      // Full lY gradient through R:
      //   ∂LL_f2u/∂lY(r,c) = -(invRER·R·sumYY)(Ri,Rj) + 2·(S_RER·R·e_sub)(Ri,Rj)
      // (the first term comes from u_i = R·y_colsR changing; the second from RER changing)
      auto add_lY_contrib = [&](arma::uword lyr, arma::uword lyc, double& g) {
        const arma::uword li = lyc * M.lY.n_rows + lyr;
        if (li < lY_sz && lY2R[li] != (arma::uword)-1) {
          const arma::uword Rlin = lY2R[li];
          const arma::uword Ri = Rlin % nU, Rj = Rlin / nU;
          g += invRER_R_YY(Ri, Rj) - 2.0 * SRR_e(Ri, Rj);
        }
      };

      for (std::size_t k = 0; k < p; ++k) {
        const arma::uword r   = row[k];
        const arma::uword c   = col[k];
        const bool        sym = static_cast<bool>(symmetric[k]);

        switch (block[k]) {
          case 1: {  // lY
            add_lY_contrib(r, c, grad[k]);
            if (sym && r != c) add_lY_contrib(c, r, grad[k]);
            break;
          }
          case 5: {  // e — only colsR diagonal entries contribute here
            for (arma::uword rr = 0; rr < nColsR; ++rr) {
              if (M.colsR_uvec(rr) == r) { grad[k] += RT_S_R(rr, rr); break; }
            }
            break;
          }
          default: break;
        }
      }
    }
  }

  // -----------------------------------------------------------------------
  // Step 5a: tY (blk 3) and a (blk 8)
  // -----------------------------------------------------------------------
  {
    // a: ∂Ey_i/∂a(r,0) = Binv_i[:,r]; constant Binv_c for kOmegaEta == 0
    arma::vec Binv_Ey;
    if (kOmega0) Binv_Ey = M.Binv_c.t() * s_Ey_acc;    // numEtas-vec

    // tY at colsR positions has two contributions:
    //   f2u: ∂LL_f2u/∂tY(colsR[jj]) = (R' · invRER · R · sumY_colsR)(jj)
    //     (from ∂u_i/∂tY = -R[:,jj], linear in u_i)
    //   f3:  ∂LL_f3/∂tY(colsR[jj])  = -(L2R_cache' · s_Ey_acc)(jj)
    //     (from Ey_i += L2R_cache · y_colsR_i, shifted by -tY)
    arma::vec RT_invRER_R_sy, L2R_Ey;
    if (M.hasR) {
      RT_invRER_R_sy = M.R.t() * M.invRER * M.R * sumY_colsR_vec;
      L2R_Ey         = M.L2R_cache.t() * s_Ey_acc;
    }

    for (std::size_t k = 0; k < p; ++k) {
      const arma::uword r = row[k];
      switch (block[k]) {
        case 3: {  // tY
          // A given tY can appear in BOTH nonNormal_uvec AND colsR_uvec
          // (e.g. the scaling indicator is in f3 directly and also in u=R·y_colsR).
          // Accumulate all matching contributions without short-circuiting.
          //
          // direct f3: dv_i(jj) = yi(jj) - Ey_i(jj), ∂yi(jj)/∂tY(r,0) = -1
          // ∂f3/∂dv * ∂dv/∂tY = (-s_Ey_i(jj)) * (-1) = +s_Ey_i(jj)
          for (int jj = 0; jj < pY_f3; ++jj) {
            if (M.nonNormal_uvec(jj) == r) grad[k] += s_Ey_acc(jj);
          }
          // colsR path: f2u linear + f3 via L2R_cache · y_colsR_i
          if (M.hasR) {
            for (arma::uword jj = 0; jj < nColsR; ++jj) {
              if (M.colsR_uvec(jj) == r) {
                grad[k] += RT_invRER_R_sy(jj) - L2R_Ey(jj);
              }
            }
          }
          break;
        }
        case 8: {  // a: ∂Ey_i/∂a(r,0) = Binv_i[:,r]
          if (kOmega0) grad[k] += Binv_Ey(r);
          // kOmegaEta != 0 needs per-obs Binv_i accumulation — deferred to later step
          break;
        }
        default: break;
      }
    }
  }

  // -----------------------------------------------------------------------
  // Step 5b: lX (blk 0), d (blk 4), phi (blk 16) — f3 contributions via Sigma1
  // -----------------------------------------------------------------------
  // SE_i += BG2O_i · Sigma1 · BG2O_i'.
  // Differentiating Sigma1 = phi - L1 · lX · phi w.r.t. LXPLX:
  //   dSigma1 = L1 · d(LXPLX) · L1'  (L1 = phi·lX'·invLXPLX)
  // → effective score w.r.t. LXPLX from f3:
  //   S_LXPLX_f3 = L1' · S_SE_BG2O_acc · L1
  // Then the same Step-3 formulas apply with S_LXPLX_f3 in place of S_LXPLX.
  {
    const arma::mat S_LXPLX_f3  = M.L1.t() * S_SE_BG2O_acc * M.L1; // pX × pX
    const arma::mat SlX_phi_f3  = S_LXPLX_f3 * M.lX * M.phi;
    const arma::mat lXT_S_lX_f3 = M.lX.t() * S_LXPLX_f3 * M.lX;

    for (std::size_t k = 0; k < p; ++k) {
      const arma::uword r   = row[k];
      const arma::uword c   = col[k];
      const bool        sym = static_cast<bool>(symmetric[k]);
      const double      fac = (sym && r != c) ? 2.0 : 1.0;

      switch (block[k]) {
        case 0:  grad[k] += 2.0 * SlX_phi_f3(r, c);   break;  // lX
        case 4:  grad[k] += S_LXPLX_f3(r, c);          break;  // d
        case 16: grad[k] += fac * lXT_S_lX_f3(r, c);  break;  // phi
        default: break;
      }
    }
  }

  // -----------------------------------------------------------------------
  // Step 5c: Psi (blk 7) and e non-colsR (blk 5) via Sig2_i
  // -----------------------------------------------------------------------
  // Sig2_i = Binv_i · Psi · Binv_i' + fullSig2TE
  //
  // Psi (kOmegaEta == 0):
  //   constant Binv_c → score = (Binv_c' · S_SE_acc · Binv_c)(r,c)
  //
  // e non-colsR has two paths into fullSig2TE:
  //   TE2: non-latent etas — fullSig2TE(p,p) = e(q,q) → grad += S_SE_acc(p,p)
  //   TE1: latent etas — fullSig2TE via Sigma2TE_sub and L2R_cache
  //     Sigma2TE_sub = subTE1 − diag(te1²) · BtRB  (BtRB = Beta' · invRER · Beta)
  //     L2R_cache[latentEta,:] = −subTE1 · Beta' · invRER · R
  //     ∂LL_f3/∂te1(j) = S_sub(j,j) − 2·te1(j)·(S_sub·BtRB)(j,j)
  //                       − s_Ey_lat(j)·(Beta'·invRER·R·sumY_colsR)(j)
  {
    arma::mat Binvt_S_Binv;
    if (kOmega0) Binvt_S_Binv = M.Binv_c.t() * S_SE_acc * M.Binv_c;

    arma::vec te1_grad;
    if (M.hasR && M.selectTE1_diag.n_elem > 0) {
      const arma::uword ns    = M.selectTE1_diag.n_elem;
      const arma::mat S_sub   = S_SE_acc.submat(M.latentEta_uvec, M.latentEta_uvec);
      const arma::mat BtRB    = M.Beta.t() * M.invRER * M.Beta;
      const arma::vec S_BtRB_diag = arma::diagvec(S_sub * BtRB);
      const arma::vec s_Ey_lat    = s_Ey_acc(M.latentEta_uvec);
      const arma::vec BtR_sumYR   = M.Beta.t() * M.invRER * M.R * sumY_colsR_vec;
      te1_grad.set_size(ns);
      for (arma::uword j = 0; j < ns; ++j) {
        te1_grad(j) = S_sub(j, j)
                    - 2.0 * M.subTE1(j, j) * S_BtRB_diag(j)
                    - s_Ey_lat(j) * BtR_sumYR(j);
      }
    }

    for (std::size_t k = 0; k < p; ++k) {
      const arma::uword r   = row[k];
      const arma::uword c   = col[k];
      const bool        sym = static_cast<bool>(symmetric[k]);
      const double      fac = (sym && r != c) ? 2.0 : 1.0;

      switch (block[k]) {
        case 7: {  // Psi
          if (kOmega0) grad[k] += fac * Binvt_S_Binv(r, c);
          // kOmegaEta != 0: needs per-obs Binv_i accumulation — deferred
          break;
        }
        case 5: {  // e non-colsR; colsR entries are handled in Step 4
          if (M.hasTE2) {
            for (arma::uword kk = 0; kk < M.selectTE2_diag.n_elem; ++kk) {
              if (M.selectTE2_diag(kk) == r) {
                grad[k] += S_SE_acc(M.nonLatentEta_uvec(kk), M.nonLatentEta_uvec(kk));
                break;
              }
            }
          }
          if (M.hasR && te1_grad.n_elem > 0) {
            for (arma::uword j = 0; j < M.selectTE1_diag.n_elem; ++j) {
              if (M.selectTE1_diag(j) == r) { grad[k] += te1_grad(j); break; }
            }
          }
          break;
        }
        default: break;
      }
    }
  }

  // TODO Step 5d–f: Gx, beta0, tX, Oxx, Ge, Oex

  return grad;
}


// Temporary export for testing Step 3 — will be removed once all steps are done.
// [[Rcpp::export]]
arma::vec analyticalGradQmlCpp(const Rcpp::List& submodel,
                                const arma::uvec& block,
                                const arma::uvec& row,
                                const arma::uvec& col,
                                const arma::uvec& symmetric) {
  try {
    const QMLModel  M(submodel);
    const arma::mat data_full =
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);
    return analyticalGradQmlCore(M, data_full, block, row, col, symmetric);
  } catch (const std::exception& e) {
    Rcpp::warning("analyticalGradQmlCpp exception: %s", e.what());
    arma::vec nan_grad(block.n_elem);
    nan_grad.fill(arma::datum::nan);
    return nan_grad;
  } catch (...) {
    Rcpp::warning("analyticalGradQmlCpp: unknown exception");
    arma::vec nan_grad(block.n_elem);
    nan_grad.fill(arma::datum::nan);
    return nan_grad;
  }
}


// [[Rcpp::export]]
arma::vec gradLogLikQmlCpp(const Rcpp::List& submodel,
                            const arma::uvec& block,
                            const arma::uvec& row,
                            const arma::uvec& col,
                            const arma::uvec& symmetric,
                            const double eps    = 1e-6,
                            const int ncores    = 1L) {
  try {
    ThreadSetter ts(ncores);
    QMLModel M(submodel);
    const arma::mat data_full =
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);

    auto qml_ll = [&data_full](QMLModel& Mc) -> double {
      return logLikQmlFromModel(Mc, data_full, 1);
    };

    return gradientFDQmlCore(M, qml_ll, block, row, col, symmetric, eps, ncores);
  } catch (...) {
    // Return NaN gradient so the optimizer steps away from this bad point.
    arma::vec nan_grad(block.n_elem);
    nan_grad.fill(arma::datum::nan);
    return nan_grad;
  }
}


// =============================================================
// Hybrid Hessian: FD of C++ gradient
// =============================================================

// [[Rcpp::export]]
Rcpp::List hessLogLikQmlCpp(const Rcpp::List& submodel,
                              const arma::uvec& block,
                              const arma::uvec& row,
                              const arma::uvec& col,
                              const arma::uvec& symmetric,
                              const double relStep  = 1e-4,
                              const double minAbs   = 1.0,
                              const int ncores      = 1L) {
  ThreadSetter ts(ncores);
  QMLModel M(submodel);
  const arma::mat data_full =
    Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);

  auto qml_ll = [&data_full](QMLModel& Mc) -> double {
    return logLikQmlFromModel(Mc, data_full, 1);
  };

  const std::size_t p = block.n_elem;
  const arma::vec base = getParamsQml(M, block, row, col);
  const arma::vec incr =
    arma::max(arma::abs(base), arma::vec(p).fill(minAbs)) * relStep;

  const double    f0    = logLikQmlFromModel(M, data_full, 1);
  const arma::vec grad0 = gradientFDQmlCore(M, qml_ll, block, row, col,
                                             symmetric, relStep, 1);

  arma::mat Hess(p, p, arma::fill::zeros);

  #pragma omp parallel for if(ncores > 1) schedule(static) \
    shared(M, data_full, block, row, col, symmetric, p, base, incr, grad0, Hess) \
    firstprivate(qml_ll)
  for (std::size_t j = 0; j < p; ++j) {
    try {
      QMLModel Mc = M.thread_clone();
      arma::vec pars = base;
      pars[j] += incr[j];
      setParamsQml(Mc, block, row, col, symmetric, pars);
      Mc.update_cache();
      const arma::vec grad_j = gradientFDQmlCore(Mc, qml_ll, block, row, col,
                                                  symmetric, relStep, 1);
      Hess.col(j) = (grad_j - grad0) / incr[j];
    } catch (...) {
      Hess.col(j).fill(arma::datum::nan);
    }
  }

  Hess = 0.5 * (Hess + Hess.t());

  return Rcpp::List::create(
    Rcpp::Named("mean")     = f0,
    Rcpp::Named("gradient") = grad0,
    Rcpp::Named("Hessian")  = Hess
  );
}
