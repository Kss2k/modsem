// TMB template for the modsem Laplace estimator.
//
// Random effects: Zeta[i,] = (xi_raw, zeta_eta) for observation i, drawn from
// N(0, Sigma_joint) where Sigma_joint = [[phi, covZetaXi'], [covZetaXi, psi]].
//
// DESIGN NOTES
//  - R headers define 'eval' as a macro, so .eval() must not appear in source.
//    Lazy Eigen expressions are forced to evaluate by assigning them into a
//    named matrix<Type> variable (Eigen evaluates on assignment).
//  - DATA_MATRIX / DATA_VECTOR return implementation-dependent wrapper types
//    that differ across TMB versions.  All helper functions therefore use free
//    template parameters (Skel, Idx) for the data-side arguments so that the
//    compiler deduces the exact type instead of requiring an implicit conversion.

#include <TMB.hpp>

// ----- helper: fill a rectangular matrix -----------------------------------
// Skel : type returned by DATA_MATRIX (rows x cols, scalar ≈ double)
// Idx  : type returned by DATA_IVECTOR (col-major, 0-based; -1 = fixed)
template<class Type, class Skel, class Idx>
matrix<Type> fill_mat(const Skel& skel, const Idx& idx,
                      const vector<Type>& theta) {
  int nr = skel.rows(), nc = skel.cols();
  matrix<Type> out(nr, nc);
  for (int i = 0; i < nr; i++)
    for (int j = 0; j < nc; j++)
      out(i, j) = Type(skel(i, j));
  for (int j = 0; j < nc; j++)
    for (int i = 0; i < nr; i++) {
      int pos = j * nr + i;
      if (idx(pos) >= 0) out(i, j) = theta(idx(pos));
    }
  return out;
}

// ----- helper: fill a column matrix (n x 1) from a 1-D data skeleton ------
template<class Type, class Skel, class Idx>
matrix<Type> fill_col(const Skel& skel, const Idx& idx,
                      const vector<Type>& theta) {
  int n = skel.size();
  matrix<Type> out(n, 1);
  for (int i = 0; i < n; i++) out(i, 0) = Type(skel(i));
  for (int i = 0; i < (int)idx.size(); i++)
    if (idx(i) >= 0) out(i, 0) = theta(idx(i));
  return out;
}

// ----- helper: fill a lower-triangular matrix (upper stays 0) --------------
template<class Type, class Skel, class Idx>
matrix<Type> fill_ltri(const Skel& skel, const Idx& idx,
                       const vector<Type>& theta) {
  int n = skel.rows();
  matrix<Type> out(n, n);
  out.setZero();
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++) {
      int pos = j * n + i;
      Type val = Type(skel(i, j));
      if (idx(pos) >= 0) val = theta(idx(pos));
      out(i, j) = val;
    }
  return out;
}

// ----- helper: fill a symmetric matrix (lower-tri determines upper-tri) ----
template<class Type, class Skel, class Idx>
matrix<Type> fill_sym(const Skel& skel, const Idx& idx,
                      const vector<Type>& theta) {
  int n = skel.rows();
  matrix<Type> out(n, n);
  out.setZero();
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++) {
      int pos = j * n + i;
      Type val = Type(skel(i, j));
      if (idx(pos) >= 0) val = theta(idx(pos));
      out(i, j) = val;
      if (i != j) out(j, i) = val;
    }
  return out;
}

// ----- helper: materialise a transpose into a new named matrix -------------
// Avoids lazy Transpose<> expressions in subsequent arithmetic.
template<class Type>
matrix<Type> tpose(const matrix<Type>& A) {
  matrix<Type> out(A.cols(), A.rows());
  out = A.transpose();
  return out;
}

// ----- helper: kron(I_k, xi_col), xi_col is (m x 1) -> (k*m x k) ---------
template<class Type>
matrix<Type> kron_I_col(int k, const matrix<Type>& xi_col) {
  int m = xi_col.rows();
  matrix<Type> out(k * m, k);
  out.setZero();
  for (int bi = 0; bi < k; bi++)
    for (int j  = 0; j  < m; j++)
      out(bi * m + j, bi) = xi_col(j, 0);
  return out;
}


// ===========================================================================
template<class Type>
Type objective_function<Type>::operator()() {

  DATA_INTEGER(nXi);
  DATA_INTEGER(nEta);
  DATA_INTEGER(nInd);

  DATA_MATRIX(Y);           // n x nInd  (NaN where missing)
  DATA_IMATRIX(obs_idx);    // n x max_nObs  (0-based; -1 = padding)
  DATA_IVECTOR(n_obs);      // n-vector

  DATA_IVECTOR(is_ordinal); // nInd  (1 = ordinal, 0 = continuous)
  DATA_IVECTOR(thresh_n);   // nInd  (# threshold entries incl. ±Inf boundary)
  DATA_MATRIX(thresh_mat);  // nInd x max_thresh

  DATA_MATRIX(Lambda_f);
  DATA_VECTOR(Theta_diag_f);
  DATA_VECTOR(tau_f);
  DATA_MATRIX(GammaXi_f);
  DATA_MATRIX(GammaEta_f);
  DATA_MATRIX(OmegaXiXi_f);    // (nEta*nXi) x nXi
  DATA_MATRIX(OmegaEtaXi_f);   // (nEta*nXi) x nEta
  DATA_VECTOR(alpha_f);
  DATA_VECTOR(beta0_f);
  DATA_MATRIX(A_lower_f);      // nXi x nXi lower Chol of phi
  DATA_MATRIX(psi_f);          // nEta x nEta symmetric
  DATA_MATRIX(covZetaXi_f);    // nEta x nXi

  DATA_IVECTOR(Lambda_idx);
  DATA_IVECTOR(Theta_diag_idx);
  DATA_IVECTOR(tau_idx);
  DATA_IVECTOR(GammaXi_idx);
  DATA_IVECTOR(GammaEta_idx);
  DATA_IVECTOR(OmegaXiXi_idx);
  DATA_IVECTOR(OmegaEtaXi_idx);
  DATA_IVECTOR(alpha_idx);
  DATA_IVECTOR(beta0_idx);
  DATA_IVECTOR(A_lower_idx);
  DATA_IVECTOR(psi_idx);
  DATA_IVECTOR(covZetaXi_idx);

  PARAMETER_VECTOR(theta);
  PARAMETER_MATRIX(Zeta);    // n x nLatent  random effects [xi_raw | zeta_eta]

  int n       = Y.rows();
  int nLatent = nXi + nEta;

  // reconstruct model matrices from theta
  matrix<Type> Lambda     = fill_mat(Lambda_f,      Lambda_idx,      theta);
  matrix<Type> GammaXi    = fill_mat(GammaXi_f,    GammaXi_idx,    theta);
  matrix<Type> GammaEta   = fill_mat(GammaEta_f,   GammaEta_idx,   theta);
  matrix<Type> OmegaXiXi  = fill_mat(OmegaXiXi_f,  OmegaXiXi_idx,  theta);
  matrix<Type> OmegaEtaXi = fill_mat(OmegaEtaXi_f, OmegaEtaXi_idx, theta);
  matrix<Type> cZX        = fill_mat(covZetaXi_f,  covZetaXi_idx,  theta);

  matrix<Type> beta0        = fill_col(beta0_f,      beta0_idx,      theta);
  matrix<Type> alpha        = fill_col(alpha_f,      alpha_idx,      theta);
  matrix<Type> tau          = fill_col(tau_f,        tau_idx,        theta);
  matrix<Type> Theta_diag_m = fill_col(Theta_diag_f, Theta_diag_idx, theta);

  matrix<Type> A_lower = fill_ltri(A_lower_f, A_lower_idx, theta);
  matrix<Type> phi(nXi, nXi);
  phi = A_lower * tpose(A_lower);           // forced by named-variable assignment

  matrix<Type> psi = fill_sym(psi_f, psi_idx, theta);

  // joint prior covariance  [[phi, cZX'], [cZX, psi]]
  matrix<Type> Sigma(nLatent, nLatent);
  Sigma.setZero();
  Sigma.topLeftCorner(nXi,  nXi)    = phi;
  Sigma.topRightCorner(nXi, nEta)   = tpose(cZX);
  Sigma.bottomLeftCorner(nEta,  nXi)  = cZX;
  Sigma.bottomRightCorner(nEta, nEta) = psi;

  using namespace density;
  MVNORM_t<Type> prior_nll(Sigma);          // – log p(Zeta_i | Sigma)

  matrix<Type> Ie = matrix<Type>::Identity(nEta, nEta);

  Type nll = Type(0);

  for (int i = 0; i < n; i++) {

    // copy row i of Zeta into a tmbutils::vector<Type> for MVNORM_t
    vector<Type> zeta_vec(nLatent);
    for (int k = 0; k < nLatent; k++) zeta_vec(k) = Zeta(i, k);

    nll += prior_nll(zeta_vec);

    // latent variables as column matrices
    matrix<Type> xi(nXi, 1), zeta_eta(nEta, 1);
    for (int j = 0; j < nXi;  j++) xi(j,       0) = Zeta(i, j)       + beta0(j, 0);
    for (int j = 0; j < nEta; j++) zeta_eta(j,  0) = Zeta(i, nXi + j);

    // structural model:  B * eta = rhs
    matrix<Type> KronXi  = kron_I_col(nEta, xi);     // (nEta*nXi) x nEta
    matrix<Type> KronXit = tpose(KronXi);             //  nEta x (nEta*nXi)

    matrix<Type> KX_OEX(nEta, nEta);   KX_OEX  = KronXit * OmegaEtaXi;
    matrix<Type> B(nEta, nEta);        B       = Ie - GammaEta - KX_OEX;

    matrix<Type> KX_OXX(nEta, nXi);   KX_OXX  = KronXit * OmegaXiXi;
    matrix<Type> KX_OXX_xi(nEta, 1);  KX_OXX_xi = KX_OXX * xi;
    matrix<Type> GXi_xi(nEta, 1);     GXi_xi  = GammaXi * xi;

    matrix<Type> rhs(nEta, 1);
    rhs = alpha + GXi_xi + KX_OXX_xi + zeta_eta;

    matrix<Type> eta(nEta, 1);
    eta = B.inverse() * rhs;

    // yhat = tau + Lambda * [xi; eta]
    matrix<Type> xieta(nLatent, 1);
    for (int j = 0; j < nXi;  j++) xieta(j,       0) = xi(j,  0);
    for (int j = 0; j < nEta; j++) xieta(nXi + j,  0) = eta(j, 0);

    matrix<Type> yhat_m(nInd, 1);
    yhat_m = tau + Lambda * xieta;

    // measurement likelihood
    for (int jj = 0; jj < n_obs(i); jj++) {
      int j = obs_idx(i, jj);
      if (j < 0 || j >= nInd) continue;

      Type yij  = Y(i, j);
      Type yhat = yhat_m(j, 0);
      Type sd_j = sqrt(Theta_diag_m(j, 0));

      if (is_ordinal(j) == 0) {
        nll -= dnorm(yij, yhat, sd_j, true);
      } else {
        int cat = CppAD::Integer(yij);
        int nt  = thresh_n(j);
        if (cat < 1 || cat >= nt) { nll += Type(1e10); continue; }
        Type lo = (Type(thresh_mat(j, cat - 1)) - yhat) / sd_j;
        Type hi = (Type(thresh_mat(j, cat))     - yhat) / sd_j;
        nll -= log(pnorm(hi) - pnorm(lo) + Type(1e-300));
      }
    }
  }

  REPORT(Zeta);
  return nll;
}
