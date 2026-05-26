#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "qml.h"
#include "mvnorm.h"

inline arma::mat qmlOexJacobianTerm(const arma::mat& omegaEtaXi,
                                    const arma::vec& etaMean,
                                    const int numXi,
                                    const int numEta) {
  const arma::vec oexEta = omegaEtaXi * etaMean;
  arma::mat out(numEta, numXi, arma::fill::none);

  for (int eta = 0; eta < numEta; ++eta) {
    const int offset = eta * numXi;
    for (int xi = 0; xi < numXi; ++xi)
      out(eta, xi) = oexEta(offset + xi);
  }

  return out;
}

inline arma::mat qmlJacobianEtaXi(const arma::mat& Binv,
                                  const arma::mat& gammaXi,
                                  const arma::mat& kronXi,
                                  const arma::mat& omegaXiXi2T,
                                  const arma::mat& omegaEtaXi,
                                  const arma::vec& etaMean,
                                  const bool includeOexTerm) {
  arma::mat jac = Binv * gammaXi + kronXi * omegaXiXi2T;

  if (includeOexTerm)
    jac += Binv * qmlOexJacobianTerm(omegaEtaXi, etaMean,
                                     gammaXi.n_cols, gammaXi.n_rows);

  return jac;
}

inline arma::vec qmlDmCovFromBG(const arma::mat& scoreBg,
                                const arma::mat& omegaXiXi2T,
                                const int numEta,
                                const int numXi) {
  arma::vec out(numXi, arma::fill::zeros);

  for (int xi = 0; xi < numXi; ++xi) {
    double acc = 0.0;
    for (int eta = 0; eta < numEta; ++eta) {
      const int row = eta * numXi + xi;
      for (int col = 0; col < numXi; ++col)
        acc += scoreBg(eta, col) * omegaXiXi2T(row, col);
    }
    out(xi) = acc;
  }

  return out;
}

inline double qmlVarZDerivative(const arma::mat& omega,
                                const arma::mat& sigma1,
                                const int row,
                                const int col) {
  const int ds = static_cast<int>(sigma1.n_rows);
  double out = 0.0;

  for (int k = 0; k < ds; ++k) {
    for (int s = 0; s < ds; ++s) {
      out += omega(k, s) *
        (sigma1(row, col) * sigma1(k, s) +
         sigma1(row, k)   * sigma1(col, s) +
         sigma1(row, s)   * sigma1(col, k));
    }
  }

  for (int i = 0; i < ds; ++i) {
    for (int j = 0; j < ds; ++j) {
      out += omega(i, j) *
        (sigma1(i, j)   * sigma1(row, col) +
         sigma1(i, row) * sigma1(j, col) +
         sigma1(i, col) * sigma1(j, row));
    }
  }

  return out - 2.0 * arma::trace(omega * sigma1) * sigma1(col, row);
}


inline double qmlVarZSigmaDerivative(const arma::mat& omega,
                                     const arma::mat& sigma1,
                                     const int row,
                                     const int col) {
  const int ds = static_cast<int>(sigma1.n_rows);
  double out = 0.0;
  double omegaSigmaElem = 0.0;

  for (int i = 0; i < ds; ++i)
    for (int j = 0; j < ds; ++j)
      omegaSigmaElem += omega(i, j) * sigma1(i, j);

  out += 2.0 * omega(row, col) * omegaSigmaElem;

  for (int j = 0; j < ds; ++j)
    for (int s = 0; s < ds; ++s)
      out += omega(row, j) * omega(col, s) * sigma1(j, s);

  for (int i = 0; i < ds; ++i)
    for (int k = 0; k < ds; ++k)
      out += omega(i, row) * omega(k, col) * sigma1(i, k);

  for (int j = 0; j < ds; ++j)
    for (int k = 0; k < ds; ++k)
      out += omega(row, j) * omega(k, col) * sigma1(j, k);

  for (int i = 0; i < ds; ++i)
    for (int s = 0; s < ds; ++s)
      out += omega(i, row) * omega(col, s) * sigma1(i, s);

  return out - 2.0 * arma::trace(omega * sigma1) * omega(col, row);
}


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
  const arma::mat alpha       = m["alpha"];
  const arma::mat beta0       = m["beta0"];
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
  const arma::mat omegaEtaXi  = m["omegaEtaXi"];
  const bool includeOexTerm   = Binv.n_rows > static_cast<unsigned>(numEta);

  const arma::mat varZ = varZCpp(omegaXiXi, Sigma1, numEta);
  const arma::vec trOmegaSigma = traceOmegaSigma1(omegaXiXi * Sigma1, numEta);
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
      const arma::vec mi       = beta0 + L1 * X.row(i).t();
      const arma::vec Ey_i     = Binv_t *
        ( trOmegaSigma + alpha
          + gammaXi * mi
          + kronXi_t * omegaXiXi * mi );

      const arma::mat BinvGammaXi2Omega =
        qmlJacobianEtaXi(Binv_t, gammaXi, kronXi_t, omegaXiXi2T,
                         omegaEtaXi, Ey_i, includeOexTerm);

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
      const arma::vec mi = beta0 + L1 * X.row(i).t();
      const arma::vec Ey_i = Binv *
        ( trOmegaSigma + alpha
          + gammaXi * mi
          + kronXi_t * omegaXiXi * mi );
      const arma::mat BinvGammaXi2Omega =
        qmlJacobianEtaXi(Binv, gammaXi, kronXi_t, omegaXiXi2T,
                         omegaEtaXi, Ey_i, includeOexTerm);

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
  //   7=Psi,8=a,9=beta0,10=Gx,11=Ge,12=Oxx,13=Oex,14=W,15=T,16=phi)
  arma::mat lX, lY, tX, tY, d, e, Psi, a, beta0, Gx, Ge, Oxx, Oex, W, T, phi, Ie;

  unsigned numXis, numEtas;
  int      kOmegaEta;
  bool     hasR;   // any latent etas with multiple indicators
  bool     hasTE2; // any non-latent etas needing residual variance in fullSig2TE

  // Column indices into data.full
  arma::uvec xColIdx;
  arma::uvec yColIdx;
  arma::vec  samplingWeights;
  bool       hasSamplingWeights;

  // Structural indices (hasR)
  arma::uvec rNaLinidx;      // col-major NA positions in emptyR
  arma::uvec lYFreeLinidx;   // col-major !selectScalingY positions in lY
  arma::uvec colsRUvec;      // positions of colsR in allIndsEtas (= e diag idx = y col idx)
  arma::uvec betaRowsUvec;   // row positions of selectBetaRows in lY
  arma::uvec latentEtaUvec;  // positions of latentEtas in etas ordering
  arma::uvec selectTE1Diag; // diagonal positions in e for scaling-indicator residuals

  // Structural indices (hasTE2)
  arma::uvec selectTE2Diag;    // diagonal positions in e for non-latent scaling residuals
  arma::uvec nonLatentEtaUvec; // positions of non-latent etas in etas ordering

  // Indices for f3: nonNormalInds = allIndsEtas NOT in colsU (scaling/single indicators)
  arma::uvec nonNormalUvec;    // positions within etasCv (= y column ordering)

  // Template for R (0 where fixed, will be overwritten at NA positions)
  arma::mat rTemplate;

  // Cached obs-independent derived quantities
  arma::vec  tauXAdj;
  arma::mat  LXPLX, cholLXPLX, invLXPLX, L1, Sigma1;    // cholLXPLX is upper-triangular
  arma::mat  R, RER, cholRER, invRER;           // cholRER is upper-triangular
  arma::mat  Beta, subTE1, L2Sub, Sigma2TESub;
  arma::mat  fullSig2TE;                        // numEtas x numEtas
  arma::mat  L2RCache;                          // numEtas x nColsR
  arma::vec  trOmSig;
  arma::mat  varZ, Oxx2T;
  arma::mat  BinvConst, Sigma2Const, BinvVarZConst; // kOmegaEta == 0 only

  void updateCache() {
    tauXAdj = tX + lX * beta0;

    LXPLX     = lX * phi * lX.t() + d;
    cholLXPLX = arma::chol(LXPLX, "upper");
    invLXPLX  = arma::inv_sympd(LXPLX);
    L1        = phi * lX.t() * invLXPLX;
    Sigma1    = phi - L1 * lX * phi;

    if (hasR) {
      R = rTemplate;
      for (arma::uword k = 0; k < rNaLinidx.n_elem; ++k)
        R(rNaLinidx(k)) = -lY(lYFreeLinidx(k));

      RER    = R * e.submat(colsRUvec, colsRUvec) * R.t();
      cholRER = arma::chol(RER, "upper");
      invRER = arma::inv_sympd(RER);

      Beta = lY.submat(betaRowsUvec, latentEtaUvec);

      const arma::uword ns = selectTE1Diag.n_elem;
      arma::vec te1(ns);
      for (arma::uword k = 0; k < ns; ++k)
        te1(k) = e(selectTE1Diag(k), selectTE1Diag(k));
      subTE1 = arma::diagmat(te1);

      L2Sub       = -subTE1 * Beta.t() * invRER;
      Sigma2TESub = subTE1 - arma::diagmat(arma::square(te1)) * Beta.t() * invRER * Beta;

      fullSig2TE.zeros();
      fullSig2TE.submat(latentEtaUvec, latentEtaUvec) = Sigma2TESub;
      L2RCache.zeros();
      L2RCache.rows(latentEtaUvec) = L2Sub * R;
    }

    if (hasTE2) {
      for (arma::uword k = 0; k < selectTE2Diag.n_elem; ++k) {
        const arma::uword p = nonLatentEtaUvec(k);
        fullSig2TE(p, p) = e(selectTE2Diag(k), selectTE2Diag(k));
      }
    }

    trOmSig = traceOmegaSigma1(Oxx * Sigma1, numEtas);
    varZ    = varZCpp(Oxx, Sigma1, numEtas);
    Oxx2T   = Oxx + transposeOmega(Oxx, numEtas);

    if (kOmegaEta == 0) {
      BinvConst     = arma::inv(Ie - Ge);
      Sigma2Const   = BinvConst * Psi * BinvConst.t() + fullSig2TE;
      BinvVarZConst = BinvConst * varZ * BinvConst.t();
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
    W     = Rcpp::as<arma::mat>(mat["W"]);
    T     = Rcpp::as<arma::mat>(mat["T"]);
    phi   = Rcpp::as<arma::mat>(mat["phi"]);
    Ie    = Rcpp::as<arma::mat>(mat["Ieta"]);

    const Rcpp::CharacterVector xisCv  = info["allIndsXis"];
    const Rcpp::CharacterVector etasCv = info["allIndsEtas"];
    const Rcpp::CharacterVector etaLv  = info["etas"];  // latent variable names

    const Rcpp::CharacterVector dataCols = qmlMatColnames(dat["data.full"]);
    xColIdx = qmlCharMatch(xisCv,  dataCols);
    yColIdx = qmlCharMatch(etasCv, dataCols);

    // nonNormalUvec: positions in etasCv NOT absorbed into u (i.e. not in colsU).
    // These are the scaling / single-indicator etas used for the f3 density.
    {
      const SEXP colsUSexp = mat["colsU"];
      if (!Rf_isNull(colsUSexp) && Rf_length(colsUSexp) > 0) {
        const Rcpp::CharacterVector colsUCv =
          Rcpp::as<Rcpp::CharacterVector>(colsUSexp);
        arma::uvec isU(etasCv.size(), arma::fill::zeros);
        for (int j = 0; j < etasCv.size(); ++j)
          for (int k = 0; k < colsUCv.size(); ++k)
            if (etasCv[j] == colsUCv[k]) { isU(j) = 1; break; }
        nonNormalUvec = arma::find(isU == 0);
      } else {
        nonNormalUvec =
          arma::regspace<arma::uvec>(0, (arma::uword)etasCv.size() - 1);
      }
    }

    hasSamplingWeights = !Rf_isNull(dat["weights"]);
    if (hasSamplingWeights) samplingWeights = Rcpp::as<arma::vec>(dat["weights"]);

    fullSig2TE.zeros(numEtas, numEtas);

    const SEXP emptyRSexp = mat["emptyR"];
    hasR = !Rf_isNull(emptyRSexp);

    if (hasR) {
      const arma::mat emptyRMat = Rcpp::as<arma::mat>(emptyRSexp);
      rNaLinidx = arma::find_nonfinite(emptyRMat);
      rTemplate = emptyRMat;
      rTemplate.replace(arma::datum::nan, 0.0);

      lYFreeLinidx = arma::find(
          Rcpp::as<arma::mat>(mat["selectScalingY"]) < 0.5);

      colsRUvec =
          qmlCharMatch(Rcpp::as<Rcpp::CharacterVector>(mat["colsR"]), etasCv);

      betaRowsUvec = arma::find(
          Rcpp::as<arma::mat>(mat["selectBetaRows"]).col(0) > 0.5);

      const Rcpp::CharacterVector latEtaCv = info["latentEtas"];
      latentEtaUvec = qmlCharMatch(latEtaCv, etaLv);

      selectTE1Diag = arma::find(
          Rcpp::as<arma::mat>(mat["selectThetaEpsilon1"]).diag() > 0.5);

      L2RCache.zeros(numEtas, colsRUvec.n_elem);
    }

    hasTE2 = false;
    const SEXP selTE2Sexp = mat["selectThetaEpsilon2"];
    if (!Rf_isNull(selTE2Sexp)) {
      const arma::mat selTE2 = Rcpp::as<arma::mat>(selTE2Sexp);
      if (selTE2.n_elem > 0) {
        selectTE2Diag = arma::find(selTE2.diag() > 0.5);
        hasTE2 = (selectTE2Diag.n_elem > 0);
      }
    }

    if (hasTE2) {
      if (hasR) {
        arma::uvec isLat(numEtas, arma::fill::zeros);
        for (arma::uword k = 0; k < latentEtaUvec.n_elem; ++k)
          isLat(latentEtaUvec(k)) = 1;
        nonLatentEtaUvec = arma::find(isLat == 0);
      } else {
        nonLatentEtaUvec = arma::regspace<arma::uvec>(0, numEtas - 1);
      }
    }

    updateCache();
  }

  QMLModel threadClone() const {
    QMLModel c = *this;
    c.lX = arma::mat(lX);  c.lY = arma::mat(lY);
    c.tX = arma::mat(tX);  c.tY = arma::mat(tY);
    c.d  = arma::mat(d);   c.e  = arma::mat(e);
    c.Psi   = arma::mat(Psi);   c.a     = arma::mat(a);
    c.beta0 = arma::mat(beta0); c.Gx    = arma::mat(Gx);
    c.Ge    = arma::mat(Ge);    c.Oxx   = arma::mat(Oxx);
    c.Oex   = arma::mat(Oex);   c.W     = arma::mat(W);
    c.T     = arma::mat(T);     c.phi   = arma::mat(phi);
    return c;
  }
};


inline double& qmlParam(QMLModel& M, std::size_t blk,
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
    case 14: return M.W    (r,c);
    case 15: return M.T    (r,c);
    case 16: return M.phi  (r,c);
    default: Rcpp::stop("qmlParam: unknown block %zu", blk);
  }
}


inline double logLikQmlFromModel(const QMLModel& M,
                                 const arma::mat& dataFull,
                                 const int ncores = 1) {
  const int t     = static_cast<int>(dataFull.n_rows);
  const int pX    = static_cast<int>(M.xColIdx.n_elem);
  const int pY    = static_cast<int>(M.yColIdx.n_elem);
  const int pYF3  = static_cast<int>(M.nonNormalUvec.n_elem); // scaling/single inds
  static const double log2piD = std::log(2.0 * M_PI);

  arma::mat X(t, pX), Y(t, pY);
  for (int j = 0; j < pX; ++j)
    X.col(j) = dataFull.col(M.xColIdx(j)) - M.tauXAdj(j);
  for (int j = 0; j < pY; ++j)
    Y.col(j) = dataFull.col(M.yColIdx(j)) - M.tY(j, 0);

  // Pre-compute constant parts of f2 log-densities
  const double logdet_x  = 2.0 * arma::sum(arma::log(M.cholLXPLX.diag()));
  const double f2xConst = -0.5 * (pX * log2piD + logdet_x);

  int    nU = 0;
  double f2uConst = 0.0;
  arma::mat yColsR;
  if (M.hasR) {
    nU = static_cast<int>(M.RER.n_rows);
    const double logdet_u = 2.0 * arma::sum(arma::log(M.cholRER.diag()));
    f2uConst = -0.5 * (nU * log2piD + logdet_u);

    const arma::uword nColsR = M.colsRUvec.n_elem;
    yColsR.set_size(t, nColsR);
    for (arma::uword j = 0; j < nColsR; ++j)
      yColsR.col(j) = Y.col(M.colsRUvec(j));
  }

  const bool kOmega0 = (M.kOmegaEta == 0);
  double ll   = 0.0;
  bool   fail = false;

#ifdef _OPENMP
  const bool inParallel = omp_in_parallel();
  const int  nth = inParallel ? 1 : ncores;
#else
  const int  nth = 1;
#endif

#pragma omp parallel for if(nth > 1) schedule(static) \
  reduction(+:ll) reduction(||:fail)
  for (int i = 0; i < t; ++i) {
    const arma::vec xi = X.row(i).t();
    // yi: only scaling/single indicators (nonNormalInds), one per eta
    arma::vec yi(pYF3);
    for (int jj = 0; jj < pYF3; ++jj) yi(jj) = Y(i, M.nonNormalUvec(jj));

    // --- f2x: logN(xi | 0, LXPLX) ---
    // cholLXPLX is upper triangular U: LXPLX = U^T U
    // maha = ||U^{-T} xi||^2, solved via trimatl(U^T)
    const arma::vec zx = arma::solve(arma::trimatl(M.cholLXPLX.t()), xi,
                                     arma::solve_opts::fast);
    const double f2x_i = f2xConst - 0.5 * arma::dot(zx, zx);

    // --- f2u: logN(u_raw | 0, RER) ---
    double f2u_i = 0.0;
    if (M.hasR) {
      const arma::vec uRaw = M.R * yColsR.row(i).t();
      const arma::vec zu   = arma::solve(arma::trimatl(M.cholRER.t()), uRaw,
                                          arma::solve_opts::fast);
      f2u_i = f2uConst - 0.5 * arma::dot(zu, zu);
    }

    // --- f3: logN(yi | Ey_i, SigmaE_i) ---
    const arma::vec mi = M.beta0 + M.L1 * xi;
    const arma::mat Ki = arma::kron(M.Ie, mi.t());

    arma::mat Binv_i, Sig2_i, BvZ_i;
    if (kOmega0) {
      Binv_i = M.BinvConst;
      Sig2_i = M.Sigma2Const;
      BvZ_i  = M.BinvVarZConst;
    } else {
      Binv_i = arma::inv(M.Ie - M.Ge - Ki * M.Oex);
      Sig2_i = Binv_i * M.Psi * Binv_i.t() + M.fullSig2TE;
      BvZ_i  = Binv_i * M.varZ * Binv_i.t();
    }

    arma::vec Ey_i = Binv_i * (M.trOmSig + M.a + M.Gx * mi + Ki * M.Oxx * mi);
    if (M.hasR)
      Ey_i += M.L2RCache * yColsR.row(i).t();

    const arma::mat BG2O =
      qmlJacobianEtaXi(Binv_i, M.Gx, Ki, M.Oxx2T, M.Oex, Ey_i, !kOmega0);
    const arma::mat SE_i = BG2O * M.Sigma1 * BG2O.t() + Sig2_i + BvZ_i;

    arma::mat Lf;
    if (!arma::chol(Lf, SE_i, "lower")) { fail = true; continue; }

    const arma::vec dv   = yi - Ey_i;
    const arma::vec zf   = arma::solve(arma::trimatl(Lf), dv, arma::solve_opts::fast);
    const double logdet_f = 2.0 * arma::sum(arma::log(Lf.diag()));
    const double f3_i    = -0.5 * (pYF3 * log2piD + logdet_f + arma::dot(zf, zf));

    const double w = M.hasSamplingWeights ? M.samplingWeights(i) : 1.0;
    ll += w * (f2x_i + f2u_i + f3_i);
  }

  if (fail) return arma::datum::nan;
  return ll;
}

inline arma::vec obsLogLikQmlFromModel(const QMLModel& M,
                                       const arma::mat& dataFull,
                                       const int ncores = 1) {
  const int t    = static_cast<int>(dataFull.n_rows);
  const int pX   = static_cast<int>(M.xColIdx.n_elem);
  const int pY   = static_cast<int>(M.yColIdx.n_elem);
  const int pYF3 = static_cast<int>(M.nonNormalUvec.n_elem);
  static const double log2piD = std::log(2.0 * M_PI);

  arma::vec ll(t, arma::fill::zeros);
  arma::mat X(t, pX), Y(t, pY);
  for (int j = 0; j < pX; ++j)
    X.col(j) = dataFull.col(M.xColIdx(j)) - M.tauXAdj(j);
  for (int j = 0; j < pY; ++j)
    Y.col(j) = dataFull.col(M.yColIdx(j)) - M.tY(j, 0);

  const double logdet_x  = 2.0 * arma::sum(arma::log(M.cholLXPLX.diag()));
  const double f2xConst = -0.5 * (pX * log2piD + logdet_x);

  int    nU = 0;
  double f2uConst = 0.0;
  arma::mat yColsR;
  if (M.hasR) {
    nU = static_cast<int>(M.RER.n_rows);
    const double logdet_u = 2.0 * arma::sum(arma::log(M.cholRER.diag()));
    f2uConst = -0.5 * (nU * log2piD + logdet_u);

    const arma::uword nColsR = M.colsRUvec.n_elem;
    yColsR.set_size(t, nColsR);
    for (arma::uword j = 0; j < nColsR; ++j)
      yColsR.col(j) = Y.col(M.colsRUvec(j));
  }

  const bool kOmega0 = (M.kOmegaEta == 0);
  bool fail = false;

  for (int i = 0; i < t; ++i) {
    const arma::vec xi = X.row(i).t();

    const arma::vec zx = arma::solve(arma::trimatl(M.cholLXPLX.t()), xi,
                                     arma::solve_opts::fast);
    const double f2x_i = f2xConst - 0.5 * arma::dot(zx, zx);

    double f2u_i = 0.0;
    if (M.hasR) {
      const arma::vec uRaw = M.R * yColsR.row(i).t();
      const arma::vec zu = arma::solve(arma::trimatl(M.cholRER.t()), uRaw,
                                       arma::solve_opts::fast);
      f2u_i = f2uConst - 0.5 * arma::dot(zu, zu);
    }

    arma::vec yi(pYF3);
    for (int jj = 0; jj < pYF3; ++jj) yi(jj) = Y(i, M.nonNormalUvec(jj));

    const arma::vec mi = M.beta0 + M.L1 * xi;
    const arma::mat Ki = arma::kron(M.Ie, mi.t());

    arma::mat Binv_i, Sig2_i, BvZ_i;
    if (kOmega0) {
      Binv_i = M.BinvConst; Sig2_i = M.Sigma2Const; BvZ_i = M.BinvVarZConst;
    } else {
      Binv_i = arma::inv(M.Ie - M.Ge - Ki * M.Oex);
      Sig2_i = Binv_i * M.Psi * Binv_i.t() + M.fullSig2TE;
      BvZ_i  = Binv_i * M.varZ * Binv_i.t();
    }

    const arma::vec h_i = M.trOmSig + M.a + M.Gx * mi + Ki * M.Oxx * mi;
    arma::vec Ey_i = Binv_i * h_i;
    if (M.hasR) Ey_i += M.L2RCache * yColsR.row(i).t();

    const arma::mat BG2O =
      qmlJacobianEtaXi(Binv_i, M.Gx, Ki, M.Oxx2T, M.Oex, Ey_i, !kOmega0);
    const arma::mat SE_i = BG2O * M.Sigma1 * BG2O.t() + Sig2_i + BvZ_i;

    arma::mat Lf;
    if (!arma::chol(Lf, SE_i, "lower")) { fail = true; continue; }
    const double logdet_f = 2.0 * arma::sum(arma::log(Lf.diag()));
    const arma::vec zf = arma::solve(arma::trimatl(Lf), yi - Ey_i,
                                     arma::solve_opts::fast);
    const double f3_i = -0.5 * (pYF3 * log2piD + logdet_f + arma::dot(zf, zf));

    const double w = M.hasSamplingWeights ? M.samplingWeights(i) : 1.0;
    ll(i) = w * (f2x_i + f2u_i + f3_i);
  }

  if (fail) ll.fill(arma::datum::nan);
  return ll;
}


// [[Rcpp::export]]
double logLikQmlCpp(const Rcpp::List& submodel, const int ncores = 1) {
  try {
    const QMLModel M(submodel);
    const arma::mat dataFull =
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);
    return logLikQmlFromModel(M, dataFull, ncores);
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
    pars[k] = qmlParam(M, block[k], row[k], col[k]);
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
    qmlParam(M, block[k], row[k], col[k]) = vals[k];
    if (symmetric[k] && row[k] != col[k])
      qmlParam(M, block[k], col[k], row[k]) = vals[k];
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
      QMLModel Mc = M.threadClone();
      double& ti = qmlParam(Mc, block[k], row[k], col[k]);
      double* tj = nullptr;
      if (symmetric[k] && row[k] != col[k])
        tj = &qmlParam(Mc, block[k], col[k], row[k]);
      ti += eps;
      if (tj) *tj += eps;
      Mc.updateCache();
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
    const arma::mat& dataFull,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric)
{
  const int        t      = static_cast<int>(dataFull.n_rows);
  const int        pX     = static_cast<int>(M.xColIdx.n_elem);
  const int        pY     = static_cast<int>(M.yColIdx.n_elem);
  const int        pYF3   = static_cast<int>(M.nonNormalUvec.n_elem);
  const std::size_t p     = block.n_elem;
  const bool       kOmega0 = (M.kOmegaEta == 0);

  arma::vec grad(p, arma::fill::zeros);

  // Mean-centered indicators
  arma::mat X(t, pX), Y(t, pY);
  for (int j = 0; j < pX; ++j)
    X.col(j) = dataFull.col(M.xColIdx(j)) - M.tauXAdj(j);
  for (int j = 0; j < pY; ++j)
    Y.col(j) = dataFull.col(M.yColIdx(j)) - M.tY(j, 0);

  const arma::uword nColsR = M.hasR ? M.colsRUvec.n_elem : 0;
  arma::mat yColsR;
  if (M.hasR) {
    yColsR.set_size(t, nColsR);
    for (arma::uword j = 0; j < nColsR; ++j)
      yColsR.col(j) = Y.col(M.colsRUvec(j));
  }

  // -----------------------------------------------------------------------
  // Accumulators
  // -----------------------------------------------------------------------
  double    sumW  = 0.0;
  arma::mat sumXX(pX, pX, arma::fill::zeros);
  arma::mat sumUU, sumYY;
  arma::mat l2RAcc;
  arma::vec sumYColsRVec;
  if (M.hasR) {
    sumUU.zeros(M.RER.n_rows, M.RER.n_rows);
    sumYY.zeros(M.colsRUvec.n_elem, M.colsRUvec.n_elem);
    l2RAcc.zeros(M.latentEtaUvec.n_elem, M.colsRUvec.n_elem);
    sumYColsRVec.zeros(nColsR);
  }

  // f3 score accumulators
  arma::vec sEyAcc(pYF3,        arma::fill::zeros);
  arma::vec eyForL2Acc(pYF3,    arma::fill::zeros);
  arma::mat sSEAcc(pYF3, pYF3, arma::fill::zeros);
  arma::mat sigma1Acc(M.Sigma1.n_rows, M.Sigma1.n_cols, arma::fill::zeros);
  arma::mat gxMeanAcc(M.Gx.n_rows, M.Gx.n_cols, arma::fill::zeros);
  arma::mat gxCovAcc(M.Gx.n_rows, M.Gx.n_cols, arma::fill::zeros);
  arma::vec alphaAcc(M.a.n_rows, arma::fill::zeros);
  arma::mat psiAcc(M.Psi.n_rows, M.Psi.n_cols, arma::fill::zeros);
  arma::vec mAcc(M.Gx.n_cols, arma::fill::zeros);
  arma::mat l1Acc(M.L1.n_rows, M.L1.n_cols, arma::fill::zeros);
  arma::vec tauXAdjAcc(pX, arma::fill::zeros);
  arma::mat geAcc(M.Ge.n_rows, M.Ge.n_cols, arma::fill::zeros);
  arma::mat oxxAcc(M.Oxx.n_rows, M.Oxx.n_cols, arma::fill::zeros);
  arma::mat oexAcc(M.Oex.n_rows, M.Oex.n_cols, arma::fill::zeros);

  bool fail = false;

  for (int i = 0; i < t; ++i) {
    const arma::vec xi = X.row(i).t();
    const double    w  = M.hasSamplingWeights ? M.samplingWeights(i) : 1.0;
    sumW  += w;
    sumXX += w * (xi * xi.t());

    if (M.hasR) {
      const arma::vec yc   = yColsR.row(i).t();
      const arma::vec uRaw = M.R * yc;
      sumUU += w * (uRaw * uRaw.t());
      sumYY += w * (yc * yc.t());
      sumYColsRVec += w * yc;
    }

    // f3 intermediates — same computation as logLikQmlFromModel
    arma::vec yi(pYF3);
    for (int jj = 0; jj < pYF3; ++jj) yi(jj) = Y(i, M.nonNormalUvec(jj));

    const arma::vec mi = M.beta0 + M.L1 * xi;
    const arma::mat Ki = arma::kron(M.Ie, mi.t());

    arma::mat Binv_i, Sig2_i, BvZ_i;
    if (kOmega0) {
      Binv_i = M.BinvConst; Sig2_i = M.Sigma2Const; BvZ_i = M.BinvVarZConst;
    } else {
      Binv_i = arma::inv(M.Ie - M.Ge - Ki * M.Oex);
      Sig2_i = Binv_i * M.Psi * Binv_i.t() + M.fullSig2TE;
      BvZ_i  = Binv_i * M.varZ * Binv_i.t();
    }

    const arma::vec h_i = M.trOmSig + M.a + M.Gx * mi + Ki * M.Oxx * mi;
    const arma::vec Ey_struct_i = Binv_i * h_i;
    arma::vec Ey_i = Ey_struct_i;
    if (M.hasR) Ey_i += M.L2RCache * yColsR.row(i).t();

    const arma::mat BG2O =
      qmlJacobianEtaXi(Binv_i, M.Gx, Ki, M.Oxx2T, M.Oex, Ey_i, !kOmega0);
    const arma::mat SE_i = BG2O * M.Sigma1 * BG2O.t() + Sig2_i + BvZ_i;

    arma::mat Lf;
    if (!arma::chol(Lf, SE_i, "lower")) { fail = true; continue; }
    const arma::mat iSE = arma::eye<arma::mat>(SE_i.n_rows, SE_i.n_cols);
    const arma::mat tmp = arma::solve(arma::trimatl(Lf), iSE,
                                       arma::solve_opts::fast);
    const arma::mat invSE_i = arma::solve(arma::trimatu(Lf.t()), tmp,
                                          arma::solve_opts::fast);

    const arma::vec dv_i    = yi - Ey_i;
    const arma::mat S_SE_i  = (-0.5) * (invSE_i - invSE_i * dv_i * dv_i.t() * invSE_i);
    const arma::vec score_y_i = invSE_i * dv_i;
    sEyAcc      += w * score_y_i;
    sSEAcc      += w * S_SE_i;
    sigma1Acc   += w * (BG2O.t() * S_SE_i * BG2O);
    alphaAcc    += w * (Binv_i.t() * score_y_i);
    psiAcc      += w * (Binv_i.t() * S_SE_i * Binv_i);

    if (M.hasR)
      l2RAcc += w * (score_y_i(M.latentEtaUvec) * yColsR.row(i));

    gxMeanAcc += w * ((Binv_i.t() * score_y_i) * mi.t());
    gxCovAcc  += w * (Binv_i.t() * S_SE_i * BG2O * M.Sigma1);

    if (!kOmega0) {
      const arma::mat scoreBg = 2.0 * S_SE_i * BG2O * M.Sigma1;
      const arma::mat tOex =
        qmlOexJacobianTerm(M.Oex, Ey_i, M.Gx.n_cols, M.Gx.n_rows);
      const arma::mat tBar = Binv_i.t() * scoreBg;

      arma::vec eyBar = score_y_i;
      for (arma::uword eta = 0; eta < tBar.n_rows; ++eta) {
        const arma::uword rowOffset = eta * M.Gx.n_cols;
        for (arma::uword xiCol = 0; xiCol < tBar.n_cols; ++xiCol) {
          const arma::uword oexRow = rowOffset + xiCol;
          const double tScore = tBar(eta, xiCol);
          oexAcc.row(oexRow) += w * tScore * Ey_i.t();
          eyBar += tScore * M.Oex.row(oexRow).t();
        }
      }

      const arma::vec eyExtra = eyBar - score_y_i;
      const arma::vec hExtra = Binv_i.t() * eyExtra;
      alphaAcc += w * hExtra;
      gxMeanAcc += w * (hExtra * mi.t());
      if (M.hasR)
        l2RAcc += w * (eyExtra(M.latentEtaUvec) * yColsR.row(i));
      eyForL2Acc += w * eyBar;

      const arma::mat Q = M.Psi + M.varZ;
      arma::mat binvBar = eyBar * h_i.t();
      binvBar += 2.0 * S_SE_i * Binv_i * Q;
      binvBar += scoreBg * (M.Gx + tOex).t();

      const arma::vec hBar = Binv_i.t() * eyBar;
      const arma::mat bBar = -Binv_i.t() * binvBar * Binv_i.t();
      const arma::mat kBar =
        scoreBg * M.Oxx2T.t() - bBar * M.Oex.t() +
        hBar * (M.Oxx * mi).t();

      geAcc += w * (-bBar);
      oexAcc += w * (-Ki.t() * bBar);

      arma::vec mBar = M.Gx.t() * hBar;
      mBar += M.Oxx.t() * (Ki.t() * hBar);
      for (arma::uword eta = 0; eta < M.Gx.n_rows; ++eta) {
        const arma::uword rowOffset = eta * M.Gx.n_cols;
        for (arma::uword xiCol = 0; xiCol < M.Gx.n_cols; ++xiCol)
          mBar(xiCol) += kBar(eta, rowOffset + xiCol);
      }

      mAcc       += w * mBar;
      l1Acc      += w * (mBar * xi.t());
      tauXAdjAcc += w * (M.invLXPLX * xi - M.L1.t() * mBar);

      const arma::mat varZBar = Binv_i.t() * S_SE_i * Binv_i;
      const arma::vec qBar = Ki.t() * hBar;
      const int numXi = static_cast<int>(M.Gx.n_cols);
      const int numEta = static_cast<int>(M.Gx.n_rows);
      for (int eta = 0; eta < numEta; ++eta) {
        const int rowOffset = eta * numXi;
        const arma::mat omegaEta =
          M.Oxx.submat(rowOffset, 0, rowOffset + numXi - 1, numXi - 1);
        sigma1Acc += w * hBar(eta) * omegaEta.t();

        for (int a = 0; a < numXi; ++a) {
          for (int b = 0; b < numXi; ++b) {
            oxxAcc(rowOffset + a, b) += w * (
              hBar(eta) * M.Sigma1(b, a) +
              qBar(rowOffset + a) * mi(b) +
              scoreBg(eta, b) * mi(a) +
              scoreBg(eta, a) * mi(b) +
              varZBar(eta, eta) *
                qmlVarZDerivative(omegaEta, M.Sigma1, a, b)
            );
            sigma1Acc(a, b) += w * varZBar(eta, eta) *
              qmlVarZSigmaDerivative(omegaEta, M.Sigma1, a, b);
          }
        }
      }
    }

    if (kOmega0) {
      arma::vec sM_i = BG2O.t() * score_y_i;
      sM_i += qmlDmCovFromBG(2.0 * S_SE_i * BG2O * M.Sigma1,
                              M.Oxx2T, M.Gx.n_rows, M.Gx.n_cols);
      eyForL2Acc += w * score_y_i;

      mAcc       += w * sM_i;
      l1Acc      += w * (sM_i * xi.t());
      tauXAdjAcc += w * (M.invLXPLX * xi - M.L1.t() * sM_i);

      const arma::mat scoreBg = 2.0 * S_SE_i * BG2O * M.Sigma1;
      const arma::mat BinvGx   = Binv_i * M.Gx;
      const arma::mat BinvQBt  = Binv_i * (M.Psi + M.varZ) * Binv_i.t();
      geAcc += w * (
        (Binv_i.t() * score_y_i) * Ey_struct_i.t() +
        Binv_i.t() * scoreBg * BinvGx.t() +
        2.0 * Binv_i.t() * S_SE_i * BinvQBt
      );

      const arma::vec sH_i = Binv_i.t() * score_y_i;
      const arma::mat varZBar = Binv_i.t() * S_SE_i * Binv_i;
      const int numXi = static_cast<int>(M.Gx.n_cols);
      const int numEta = static_cast<int>(M.Gx.n_rows);
      for (int eta = 0; eta < numEta; ++eta) {
        const int rowOffset = eta * numXi;
        const arma::mat omegaEta =
          M.Oxx.submat(rowOffset, 0, rowOffset + numXi - 1, numXi - 1);
        sigma1Acc += w * sH_i(eta) * omegaEta.t();

        for (int a = 0; a < numXi; ++a) {
          for (int b = 0; b < numXi; ++b) {
            oxxAcc(rowOffset + a, b) += w * (
              sH_i(eta) * (M.Sigma1(b, a) + mi(a) * mi(b)) +
              scoreBg(eta, b) * mi(a) +
              scoreBg(eta, a) * mi(b)
            );
            sigma1Acc(a, b) += w * varZBar(eta, eta) *
              qmlVarZSigmaDerivative(omegaEta, M.Sigma1, a, b);
          }
        }
      }
    }
  }

  if (fail) { grad.fill(arma::datum::nan); return grad; }

  // -----------------------------------------------------------------------
  // f2x and f3 exogenous-side reverse pass
  // -----------------------------------------------------------------------
  // The exogenous measurement parameters enter f3 through both:
  //   Sigma1 = phi - L1 * lX * phi
  //   L1     = phi * lX' * inv(lX * phi * lX' + thetaDelta)
  // Earlier code only pushed the Sigma1 score through LXPLX, which missed
  // the posterior-mean path mi = beta0 + L1 * x.
  arma::mat sLXPLX = -0.5 * sumW * M.invLXPLX
                     + 0.5 * M.invLXPLX * sumXX * M.invLXPLX;
  arma::mat lXAcc(M.lX.n_rows, M.lX.n_cols, arma::fill::zeros);
  arma::mat phiAcc(M.phi.n_rows, M.phi.n_cols, arma::fill::zeros);

  {
    // Reverse Sigma1 = phi - L1 * lX * phi.
    phiAcc += sigma1Acc;
    l1Acc  += -sigma1Acc * (M.lX * M.phi).t();
    lXAcc  += -M.L1.t() * sigma1Acc * M.phi.t();
    phiAcc += -(M.L1 * M.lX).t() * sigma1Acc;

    // Reverse L1 = phi * lX' * invLXPLX.
    const arma::mat h = M.phi * M.lX.t();
    phiAcc += l1Acc * M.invLXPLX.t() * M.lX;
    lXAcc  += M.invLXPLX * l1Acc.t() * M.phi;
    sLXPLX += -M.invLXPLX * h.t() * l1Acc * M.invLXPLX;

    // Reverse LXPLX = lX * phi * lX' + thetaDelta.
    const arma::mat sLXPLXSym = 0.5 * (sLXPLX + sLXPLX.t());
    lXAcc  += 2.0 * sLXPLXSym * M.lX * M.phi;
    phiAcc += M.lX.t() * sLXPLXSym * M.lX;

    for (std::size_t k = 0; k < p; ++k) {
      const arma::uword r   = row[k];
      const arma::uword c   = col[k];
      const bool        sym = static_cast<bool>(symmetric[k]);

      switch (block[k]) {
        case 0: {
          grad[k] += lXAcc(r, c);
          break;
        }
        case 4: {
          grad[k] += sLXPLXSym(r, c);
          if (sym && r != c) grad[k] += sLXPLXSym(c, r);
          break;
        }
        case 16: {
          grad[k] += phiAcc(r, c);
          if (sym && r != c) grad[k] += phiAcc(c, r);
          break;
        }
        default: break;
      }
    }
  }

  // -----------------------------------------------------------------------
  // f2u: score w.r.t. RER
  // -----------------------------------------------------------------------
  if (M.hasR) {
    const arma::mat sRER = -0.5 * sumW * M.invRER
                          + 0.5 * M.invRER * sumUU * M.invRER;
    {
      const arma::mat eSub = M.e.submat(M.colsRUvec, M.colsRUvec);
      const arma::uword nColsR = M.colsRUvec.n_elem;
      const arma::uword nU     = M.R.n_rows;

      arma::mat gR = -M.invRER * M.R * sumYY
                    + 2.0 * sRER * M.R * eSub;
      arma::mat gESub = M.R.t() * sRER * M.R;

      if (M.selectTE1Diag.n_elem > 0) {
        const arma::mat gSigma =
          sSEAcc.submat(M.latentEtaUvec, M.latentEtaUvec);
        const arma::mat D = M.subTE1;
        const arma::mat D2 = arma::diagmat(arma::square(D.diag()));
        const arma::mat C = M.invRER;
        const arma::mat A = M.Beta;
        const arma::mat H = A.t() * C * A;
        const arma::mat L2 = -D * A.t() * C;

        arma::mat gL2 = l2RAcc * M.R.t();
        gR += L2.t() * l2RAcc;

        arma::vec gTe(D.n_rows, arma::fill::zeros);
        const arma::mat AtC = A.t() * C;
        for (arma::uword j = 0; j < D.n_rows; ++j) {
          gTe(j) += gSigma(j, j);
          gTe(j) -= 2.0 * D(j, j) * arma::dot(gSigma.row(j), H.row(j));
          gTe(j) -= arma::dot(gL2.row(j), AtC.row(j));
        }

        arma::mat gH = -D2 * gSigma;
        arma::mat gA = C * A * gH.t() + C.t() * A * gH;
        arma::mat gC = A * gH.t() * A.t();

        gA += -(C * gL2.t() * D);
        gC += -(gL2.t() * D * A.t());

        const arma::mat gRER = -C.t() * gC * C.t();
        gR     += gRER * M.R * eSub + gRER.t() * M.R * eSub;
        gESub += M.R.t() * gRER * M.R;

        for (arma::uword j = 0; j < M.selectTE1Diag.n_elem; ++j) {
          const arma::uword idx = M.selectTE1Diag(j);
          gESub(idx, idx) += gTe(j);
        }

        for (std::size_t k = 0; k < p; ++k) {
          if (block[k] != 1) continue;
          const arma::uword r = row[k], c = col[k];
          const bool sym = static_cast<bool>(symmetric[k]);
          if (r < M.lY.n_rows && c < M.lY.n_cols) {
            for (arma::uword br = 0; br < M.betaRowsUvec.n_elem; ++br)
              for (arma::uword bc = 0; bc < M.latentEtaUvec.n_elem; ++bc)
                if (M.betaRowsUvec(br) == r && M.latentEtaUvec(bc) == c)
                  grad[k] += gA(br, bc);
          }
          if (sym && r != c && c < M.lY.n_rows && r < M.lY.n_cols) {
            for (arma::uword br = 0; br < M.betaRowsUvec.n_elem; ++br)
              for (arma::uword bc = 0; bc < M.latentEtaUvec.n_elem; ++bc)
                if (M.betaRowsUvec(br) == c && M.latentEtaUvec(bc) == r)
                  grad[k] += gA(br, bc);
          }
        }
      }

      // Lookup: lY linear col-major index → R linear col-major index.
      // Unused entries stay at sentinel (arma::uword)-1.
      const arma::uword lYSize = M.lY.n_rows * M.lY.n_cols;
      std::vector<arma::uword> lY2R(lYSize, (arma::uword)-1);
      for (arma::uword kk = 0; kk < M.lYFreeLinidx.n_elem; ++kk)
        lY2R[M.lYFreeLinidx(kk)] = M.rNaLinidx(kk);

      auto addLYContrib = [&](arma::uword lyr, arma::uword lyc, double& g) {
        const arma::uword li = lyc * M.lY.n_rows + lyr;
        if (li < lYSize && lY2R[li] != (arma::uword)-1) {
          const arma::uword Rlin = lY2R[li];
          const arma::uword Ri = Rlin % nU, Rj = Rlin / nU;
          g -= gR(Ri, Rj);
        }
      };

      for (std::size_t k = 0; k < p; ++k) {
        const arma::uword r   = row[k];
        const arma::uword c   = col[k];
        const bool        sym = static_cast<bool>(symmetric[k]);

        switch (block[k]) {
          case 1: {  // lY
            addLYContrib(r, c, grad[k]);
            if (sym && r != c) addLYContrib(c, r, grad[k]);
            break;
          }
          case 5: {  // e — only colsR diagonal entries contribute here
            for (arma::uword rr = 0; rr < nColsR; ++rr) {
              if (M.colsRUvec(rr) == r) { grad[k] += gESub(rr, rr); break; }
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
    // tY at colsR positions has two contributions:
    //   f2u: ∂LL_f2u/∂tY(colsR[jj]) = (R' · invRER · R · sumYColsR)(jj)
    //     (from ∂u_i/∂tY = -R[:,jj], linear in u_i)
    //   f3:  ∂LL_f3/∂tY(colsR[jj])  = -(L2RCache' · sEyAcc)(jj)
    //     (from Ey_i += L2RCache · yColsR_i, shifted by -tY)
    arma::vec rtInvRerRSy, l2REy;
    if (M.hasR) {
      rtInvRerRSy = M.R.t() * M.invRER * M.R * sumYColsRVec;
      l2REy       = M.L2RCache.t() * eyForL2Acc;
    }

    for (std::size_t k = 0; k < p; ++k) {
      const arma::uword r = row[k];
      switch (block[k]) {
        case 3: {  // tY
          // A given tY can appear in BOTH nonNormalUvec AND colsRUvec
          // (e.g. the scaling indicator is in f3 directly and also in u=R·y_colsR).
          // Accumulate all matching contributions without short-circuiting.
          //
          // direct f3: dv_i(jj) = yi(jj) - Ey_i(jj), ∂yi(jj)/∂tY(r,0) = -1
          // ∂f3/∂dv * ∂dv/∂tY = (-s_Ey_i(jj)) * (-1) = +s_Ey_i(jj)
          for (int jj = 0; jj < pYF3; ++jj) {
            if (M.nonNormalUvec(jj) == r) grad[k] += sEyAcc(jj);
          }
          // colsR path: f2u linear + f3 via L2RCache · yColsR_i
          if (M.hasR) {
            for (arma::uword jj = 0; jj < nColsR; ++jj) {
              if (M.colsRUvec(jj) == r) {
                grad[k] += rtInvRerRSy(jj) - l2REy(jj);
              }
            }
          }
          break;
        }
        case 8: {  // a: ∂Ey_i/∂a(r,0) = Binv_i[:,r]
          grad[k] += alphaAcc(r);
          break;
        }
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
  //   constant BinvConst → score = (BinvConst' · sSEAcc · BinvConst)(r,c)
  //
  // e non-colsR has one direct path into fullSig2TE:
  //   TE2: non-latent etas — fullSig2TE(p,p) = e(q,q) → grad += S_SE_acc(p,p)
  // TE1 and colsR residual paths are handled in the R/RER adjoint above.
  {
    for (std::size_t k = 0; k < p; ++k) {
      const arma::uword r   = row[k];
      const arma::uword c   = col[k];
      const bool        sym = static_cast<bool>(symmetric[k]);
      const double      fac = (sym && r != c) ? 2.0 : 1.0;

      switch (block[k]) {
        case 7: {  // Psi
          grad[k] += fac * psiAcc(r, c);
          break;
        }
        case 5: {  // e non-colsR; colsR entries are handled in Step 4
          if (M.hasTE2) {
            for (arma::uword kk = 0; kk < M.selectTE2Diag.n_elem; ++kk) {
              if (M.selectTE2Diag(kk) == r) {
                grad[k] += sSEAcc(M.nonLatentEtaUvec(kk), M.nonLatentEtaUvec(kk));
                break;
              }
            }
          }
          break;
        }
        default: break;
      }
    }
  }

  // Step 5d-f: remaining structural mappings.
  {
    const arma::vec beta0Acc = mAcc + M.lX.t() * tauXAdjAcc;
    if (kOmega0) {
      const int numXi = static_cast<int>(M.Gx.n_cols);
      const int numEta = static_cast<int>(M.Gx.n_rows);
      for (int eta = 0; eta < numEta; ++eta) {
        const arma::mat omegaEta =
          M.Oxx.submat(eta * numXi, 0, (eta + 1) * numXi - 1, numXi - 1);
        for (int a = 0; a < numXi; ++a) {
          for (int b = 0; b < numXi; ++b) {
            oxxAcc(eta * numXi + a, b) +=
              psiAcc(eta, eta) * qmlVarZDerivative(omegaEta, M.Sigma1, a, b);
          }
        }
      }
    }

    for (std::size_t k = 0; k < p; ++k) {
      const arma::uword r   = row[k];
      const arma::uword c   = col[k];
      const bool        sym = static_cast<bool>(symmetric[k]);

      switch (block[k]) {
        case 2: { // tX / tauX
          grad[k] += tauXAdjAcc(r);
          break;
        }
        case 9: { // beta0
          grad[k] += beta0Acc(r);
          break;
        }
        case 10: { // Gx / gammaXi
          grad[k] += gxMeanAcc(r, c) + 2.0 * gxCovAcc(r, c);
          if (sym && r != c)
            grad[k] += gxMeanAcc(c, r) + 2.0 * gxCovAcc(c, r);
          break;
        }
        case 11: { // Ge / gammaEta
          grad[k] += geAcc(r, c);
          if (sym && r != c) grad[k] += geAcc(c, r);
          break;
        }
        case 12: { // Oxx / omegaXiXi
          grad[k] += oxxAcc(r, c);
          if (sym && r != c && c < oxxAcc.n_rows && r < oxxAcc.n_cols)
            grad[k] += oxxAcc(c, r);
          break;
        }
        case 13: { // Oex / omegaEtaXi
          if (!kOmega0) {
            grad[k] += oexAcc(r, c);
            if (sym && r != c && c < oexAcc.n_rows && r < oexAcc.n_cols)
              grad[k] += oexAcc(c, r);
          }
          break;
        }
        default: break;
      }
    }
  }

  for (std::size_t k = 0; k < p; ++k) {
    const arma::uword blk = block[k];
    const bool supported =
      (blk == 0) ||                    // lambdaX
      (blk == 1) ||                    // lambdaY
      (blk == 2) ||                    // tauX
      (blk == 3) ||                    // tauY
      (blk == 4) ||                    // thetaDelta
      (blk == 5) ||                    // thetaEpsilon
      (blk == 9) ||                    // beta0
      (blk == 10) ||                   // gammaXi
      (blk == 11) ||                   // gammaEta
      (blk == 12) ||                   // omegaXiXi
      (blk == 13) ||                   // omegaEtaXi
      (blk == 7) ||                    // psi
      (blk == 8) ||                    // alpha
      (blk == 16);                     // phi

    if (!supported)
      grad[k] = arma::datum::nan;
  }

  return grad;
}


arma::vec hybridGradQmlCore(
    const QMLModel& M,
    const arma::mat& dataFull,
    const arma::uvec& block,
    const arma::uvec& row,
    const arma::uvec& col,
    const arma::uvec& symmetric,
    const double eps,
    const int ncores)
{
  arma::vec grad = analyticalGradQmlCore(M, dataFull, block, row, col,
                                         symmetric);
  if (grad.is_finite()) return grad;

  auto qmlLl = [&dataFull](QMLModel& Mc) -> double {
    return logLikQmlFromModel(Mc, dataFull, 1);
  };
  QMLModel Mc = M.threadClone();
  const arma::vec gradFD = gradientFDQmlCore(Mc, qmlLl, block, row, col,
                                             symmetric, eps, ncores);

  for (arma::uword k = 0; k < grad.n_elem; ++k) {
    if (!std::isfinite(grad(k))) grad(k) = gradFD(k);
  }

  return grad;
}


// [[Rcpp::export]]
arma::mat analyticalObsGradQmlCpp(const Rcpp::List& submodel,
                                  const arma::uvec& block,
                                  const arma::uvec& row,
                                  const arma::uvec& col,
                                  const arma::uvec& symmetric) {
  try {
    const QMLModel M(submodel);
    const arma::mat dataFull =
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);

    const arma::uword n = dataFull.n_rows;
    const arma::uword p = block.n_elem;
    arma::mat out(n, p, arma::fill::zeros);

    for (arma::uword i = 0; i < n; ++i) {
      QMLModel Mi = M.threadClone();
      if (Mi.hasSamplingWeights) {
        arma::vec wi(1);
        wi(0) = M.samplingWeights(i);
        Mi.samplingWeights = wi;
      }
      out.row(i) = analyticalGradQmlCore(Mi, dataFull.rows(i, i),
                                         block, row, col, symmetric).t();
    }

    return out;
  } catch (const std::exception& e) {
    Rcpp::warning("analyticalObsGradQmlCpp exception: %s", e.what());
    arma::mat nanGrad(
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]).n_rows,
      block.n_elem
    );
    nanGrad.fill(arma::datum::nan);
    return nanGrad;
  } catch (...) {
    Rcpp::warning("analyticalObsGradQmlCpp: unknown exception");
    arma::mat nanGrad(
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]).n_rows,
      block.n_elem
    );
    nanGrad.fill(arma::datum::nan);
    return nanGrad;
  }
}


// [[Rcpp::export]]
arma::vec analyticalGradQmlCpp(const Rcpp::List& submodel,
                                const arma::uvec& block,
                                const arma::uvec& row,
                                const arma::uvec& col,
                                const arma::uvec& symmetric) {
  try {
    const QMLModel  M(submodel);
    const arma::mat dataFull =
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);
    return analyticalGradQmlCore(M, dataFull, block, row, col, symmetric);
  } catch (const std::exception& e) {
    Rcpp::warning("analyticalGradQmlCpp exception: %s", e.what());
    arma::vec nanGrad(block.n_elem);
    nanGrad.fill(arma::datum::nan);
    return nanGrad;
  } catch (...) {
    Rcpp::warning("analyticalGradQmlCpp: unknown exception");
    arma::vec nanGrad(block.n_elem);
    nanGrad.fill(arma::datum::nan);
    return nanGrad;
  }
}


// [[Rcpp::export]]
arma::vec gradLogLikQmlFDCpp(const Rcpp::List& submodel,
                              const arma::uvec& block,
                              const arma::uvec& row,
                              const arma::uvec& col,
                              const arma::uvec& symmetric,
                              const double eps    = 1e-6,
                              const int ncores    = 1L) {
  try {
    ThreadSetter ts(ncores);
    QMLModel M(submodel);
    const arma::mat dataFull =
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);

    auto qmlLl = [&dataFull](QMLModel& Mc) -> double {
      return logLikQmlFromModel(Mc, dataFull, 1);
    };

    return gradientFDQmlCore(M, qmlLl, block, row, col, symmetric, eps, ncores);
  } catch (...) {
    // Return NaN gradient so the optimizer steps away from this bad point.
    arma::vec nanGrad(block.n_elem);
    nanGrad.fill(arma::datum::nan);
    return nanGrad;
  }
}


// [[Rcpp::export]]
arma::mat gradObsLogLikQmlFDCpp(const Rcpp::List& submodel,
                                const arma::uvec& block,
                                const arma::uvec& row,
                                const arma::uvec& col,
                                const arma::uvec& symmetric,
                                const double eps = 1e-6,
                                const int ncores = 1L) {
  try {
    ThreadSetter ts(ncores);
    QMLModel M(submodel);
    const arma::mat dataFull =
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);

    const arma::uword n = dataFull.n_rows;
    const arma::uword p = block.n_elem;
    arma::mat grad(n, p, arma::fill::zeros);
    const arma::vec f0 = obsLogLikQmlFromModel(M, dataFull, 1);

    #pragma omp parallel for if(ncores > 1) schedule(static) \
      shared(M, dataFull, block, row, col, symmetric, eps, grad, f0, p)
    for (std::size_t k = 0; k < p; ++k) {
      try {
        QMLModel Mc = M.threadClone();
        double& ti = qmlParam(Mc, block[k], row[k], col[k]);
        double* tj = nullptr;
        if (symmetric[k] && row[k] != col[k])
          tj = &qmlParam(Mc, block[k], col[k], row[k]);
        ti += eps;
        if (tj) *tj += eps;
        Mc.updateCache();
        grad.col(k) = (obsLogLikQmlFromModel(Mc, dataFull, 1) - f0) / eps;
      } catch (...) {
        grad.col(k).fill(arma::datum::nan);
      }
    }

    return grad;
  } catch (...) {
    arma::mat nanGrad(
      Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]).n_rows,
      block.n_elem
    );
    nanGrad.fill(arma::datum::nan);
    return nanGrad;
  }
}


// =============================================================
// Hybrid Hessian: finite difference of analytical gradients.
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
  const arma::mat dataFull =
    Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(submodel["data"])["data.full"]);

  const std::size_t p = block.n_elem;
  const arma::vec base = getParamsQml(M, block, row, col);
  const arma::vec incr =
    arma::max(arma::abs(base), arma::vec(p).fill(minAbs)) * relStep;

  const double    f0    = logLikQmlFromModel(M, dataFull, 1);
  const arma::vec grad0 = hybridGradQmlCore(
    M, dataFull, block, row, col, symmetric, relStep, 1
  );

  arma::mat Hess(p, p, arma::fill::zeros);

  #pragma omp parallel for if(ncores > 1) schedule(static) \
    shared(M, dataFull, block, row, col, symmetric, p, base, incr, grad0, Hess)
  for (std::size_t j = 0; j < p; ++j) {
    try {
      QMLModel Mc = M.threadClone();
      arma::vec pars = base;
      pars[j] += incr[j];
      setParamsQml(Mc, block, row, col, symmetric, pars);
      Mc.updateCache();
      const arma::vec gradJ = hybridGradQmlCore(
        Mc, dataFull, block, row, col, symmetric, relStep, 1
      );
      Hess.col(j) = (gradJ - grad0) / incr[j];
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
