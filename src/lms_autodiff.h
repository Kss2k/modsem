#ifndef MODSEM_LMS_AUTODIFF_H
#define MODSEM_LMS_AUTODIFF_H

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

namespace modsem {
namespace lms {

template<typename Scalar>
using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar>
using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

struct Dual2 {
  double val = 0.0;
  double d1  = 0.0;
  double d2  = 0.0;

  Dual2() = default;
  explicit Dual2(double v) : val(v), d1(0.0), d2(0.0) {}
  Dual2(double v, double d1_, double d2_) : val(v), d1(d1_), d2(d2_) {}
};

inline Dual2 make_dual2(double v, double d1 = 0.0, double d2 = 0.0) {
  return Dual2(v, d1, d2);
}

inline Dual2 operator+(const Dual2& a, const Dual2& b) {
  return Dual2(a.val + b.val, a.d1 + b.d1, a.d2 + b.d2);
}

inline Dual2 operator-(const Dual2& a, const Dual2& b) {
  return Dual2(a.val - b.val, a.d1 - b.d1, a.d2 - b.d2);
}

inline Dual2 operator-(const Dual2& a) {
  return Dual2(-a.val, -a.d1, -a.d2);
}

inline Dual2 operator*(const Dual2& a, const Dual2& b) {
  const double val = a.val * b.val;
  const double d1  = a.val * b.d1 + b.val * a.d1;
  const double d2  = a.val * b.d2 + b.val * a.d2 + 2.0 * a.d1 * b.d1;
  return Dual2(val, d1, d2);
}

inline Dual2 reciprocal(const Dual2& b) {
  const double inv_val = 1.0 / b.val;
  const double d1 = -b.d1 * inv_val * inv_val;
  const double d2 = (2.0 * b.d1 * b.d1 - b.val * b.d2) * inv_val * inv_val * inv_val;
  return Dual2(inv_val, d1, d2);
}

inline Dual2 operator/(const Dual2& a, const Dual2& b) {
  return a * reciprocal(b);
}

inline Dual2 operator+(const Dual2& a, double b) {
  return Dual2(a.val + b, a.d1, a.d2);
}
inline Dual2 operator+(double a, const Dual2& b) { return b + a; }

inline Dual2 operator-(const Dual2& a, double b) {
  return Dual2(a.val - b, a.d1, a.d2);
}
inline Dual2 operator-(double a, const Dual2& b) {
  return Dual2(a - b.val, -b.d1, -b.d2);
}

inline Dual2 operator*(const Dual2& a, double b) {
  return Dual2(a.val * b, a.d1 * b, a.d2 * b);
}
inline Dual2 operator*(double a, const Dual2& b) { return b * a; }

inline Dual2 operator/(const Dual2& a, double b) {
  const double inv = 1.0 / b;
  return Dual2(a.val * inv, a.d1 * inv, a.d2 * inv);
}

inline Dual2 operator/(double a, const Dual2& b) {
  return Dual2(a, 0.0, 0.0) / b;
}

inline Dual2& operator+=(Dual2& a, const Dual2& b) {
  a.val += b.val;
  a.d1  += b.d1;
  a.d2  += b.d2;
  return a;
}

inline Dual2& operator-=(Dual2& a, const Dual2& b) {
  a.val -= b.val;
  a.d1  -= b.d1;
  a.d2  -= b.d2;
  return a;
}

inline Dual2& operator*=(Dual2& a, double b) {
  a.val *= b;
  a.d1  *= b;
  a.d2  *= b;
  return a;
}

inline Dual2& operator*=(Dual2& a, const Dual2& b) {
  const double val = a.val * b.val;
  const double d1  = a.val * b.d1 + b.val * a.d1;
  const double d2  = a.val * b.d2 + b.val * a.d2 + 2.0 * a.d1 * b.d1;
  a.val = val;
  a.d1  = d1;
  a.d2  = d2;
  return a;
}

inline Dual2& operator/=(Dual2& a, double b) {
  const double inv = 1.0 / b;
  a.val *= inv;
  a.d1  *= inv;
  a.d2  *= inv;
  return a;
}

inline Dual2& operator/=(Dual2& a, const Dual2& b) {
  a = a / b;
  return a;
}

inline bool operator<(const Dual2& a, const Dual2& b) { return a.val < b.val; }
inline bool operator<(const Dual2& a, double b) { return a.val < b; }
inline bool operator<(double a, const Dual2& b) { return a < b.val; }

inline bool operator<=(const Dual2& a, const Dual2& b) { return a.val <= b.val; }
inline bool operator<=(const Dual2& a, double b) { return a.val <= b; }
inline bool operator<=(double a, const Dual2& b) { return a <= b.val; }

inline bool operator>(const Dual2& a, const Dual2& b) { return a.val > b.val; }
inline bool operator>(const Dual2& a, double b) { return a.val > b; }
inline bool operator>(double a, const Dual2& b) { return a > b.val; }

inline bool operator>=(const Dual2& a, const Dual2& b) { return a.val >= b.val; }
inline bool operator>=(const Dual2& a, double b) { return a.val >= b; }
inline bool operator>=(double a, const Dual2& b) { return a >= b.val; }

inline bool operator==(const Dual2& a, const Dual2& b) { return a.val == b.val; }
inline bool operator!=(const Dual2& a, const Dual2& b) { return a.val != b.val; }

inline Dual2 exp(const Dual2& x) {
  const double ex = std::exp(x.val);
  const double d1 = ex * x.d1;
  const double d2 = ex * (x.d2 + x.d1 * x.d1);
  return Dual2(ex, d1, d2);
}

inline Dual2 log(const Dual2& x) {
  const double inv = 1.0 / x.val;
  const double d1 = x.d1 * inv;
  const double d2 = (x.d2 * x.val - x.d1 * x.d1) * inv * inv;
  return Dual2(std::log(x.val), d1, d2);
}

inline Dual2 sqrt(const Dual2& x) {
  const double root = std::sqrt(x.val);
  const double d1 = 0.5 * x.d1 / root;
  const double d2 = 0.5 * (x.d2 / root) - 0.25 * (x.d1 * x.d1) / (x.val * root);
  return Dual2(root, d1, d2);
}

inline Dual2 pow(const Dual2& x, double p) {
  const double val = std::pow(x.val, p);
  const double d1 = p * std::pow(x.val, p - 1.0) * x.d1;
  const double d2 = p * std::pow(x.val, p - 1.0) * x.d2 +
    p * (p - 1.0) * std::pow(x.val, p - 2.0) * x.d1 * x.d1;
  return Dual2(val, d1, d2);
}

inline Dual2 sin(const Dual2& x) {
  const double s = std::sin(x.val);
  const double c = std::cos(x.val);
  return Dual2(s, c * x.d1, c * x.d2 - s * x.d1 * x.d1);
}

inline Dual2 cos(const Dual2& x) {
  const double s = std::sin(x.val);
  const double c = std::cos(x.val);
  return Dual2(c, -s * x.d1, -s * x.d2 - c * x.d1 * x.d1);
}

inline Dual2 tanh(const Dual2& x) {
  const double t = std::tanh(x.val);
  const double sech2 = 1.0 - t * t;
  return Dual2(t, sech2 * x.d1, sech2 * x.d2 - 2.0 * t * sech2 * x.d1 * x.d1);
}

inline Dual2 abs(const Dual2& x) {
  if (x.val >= 0.0) return x;
  return Dual2(-x.val, -x.d1, -x.d2);
}

template<typename Scalar>
struct ScalarTraits {
  static constexpr bool is_double = std::is_same_v<Scalar, double>;
  static constexpr bool is_dual2  = std::is_same_v<Scalar, Dual2>;
  static constexpr bool is_autodiff =
      !is_double && !is_dual2;
};

template<typename Scalar>
inline Eigen::Index deriv_dim_of(const Scalar& x) {
  if constexpr (ScalarTraits<Scalar>::is_autodiff) {
    return x.derivatives().size();
  } else {
    return 0;
  }
}

template<typename Scalar>
inline Eigen::Index deriv_dim_of_matrix(const Matrix<Scalar>& M) {
  if (M.size() == 0) return 0;
  return deriv_dim_of(M(0,0));
}

template<typename Scalar>
inline Scalar make_scalar(double value, Eigen::Index deriv_dim) {
  if constexpr (ScalarTraits<Scalar>::is_double) {
    return static_cast<double>(value);
  } else if constexpr (ScalarTraits<Scalar>::is_dual2) {
    return Dual2(value);
  } else {
    Scalar s(value);
    auto& d = s.derivatives();
    d.resize(deriv_dim);
    d.setZero();
    return s;
  }
}

template<typename Scalar>
inline Vector<Scalar> vec_cast(const Eigen::VectorXd& src,
                               Eigen::Index deriv_dim = 0) {
  Vector<Scalar> out(src.size());
  for (Eigen::Index i = 0; i < src.size(); ++i)
    out(i) = make_scalar<Scalar>(src(i), deriv_dim);
  return out;
}

template<typename Scalar>
inline Matrix<Scalar> mat_cast(const Eigen::MatrixXd& src,
                               Eigen::Index deriv_dim = 0) {
  Matrix<Scalar> out(src.rows(), src.cols());
  for (Eigen::Index i = 0; i < src.rows(); ++i)
    for (Eigen::Index j = 0; j < src.cols(); ++j)
      out(i, j) = make_scalar<Scalar>(src(i, j), deriv_dim);
  return out;
}

template<typename Scalar>
inline Vector<Scalar> zero_vector(Eigen::Index n, Eigen::Index deriv_dim) {
  Vector<Scalar> out(n);
  for (Eigen::Index i = 0; i < n; ++i)
    out(i) = make_scalar<Scalar>(0.0, deriv_dim);
  return out;
}

template<typename Scalar>
inline Matrix<Scalar> zero_matrix(Eigen::Index r, Eigen::Index c,
                                  Eigen::Index deriv_dim) {
  Matrix<Scalar> out(r, c);
  for (Eigen::Index i = 0; i < r; ++i)
    for (Eigen::Index j = 0; j < c; ++j)
      out(i, j) = make_scalar<Scalar>(0.0, deriv_dim);
  return out;
}

template<typename Scalar>
inline double scalar_value(const Scalar& x) {
  if constexpr (ScalarTraits<Scalar>::is_double) {
    return x;
  } else if constexpr (ScalarTraits<Scalar>::is_dual2) {
    return x.val;
  } else {
    return x.value();
  }
}

template<typename Scalar>
inline Eigen::MatrixXd to_double_matrix(const Matrix<Scalar>& src) {
  Eigen::MatrixXd out(src.rows(), src.cols());
  for (Eigen::Index i = 0; i < src.rows(); ++i)
    for (Eigen::Index j = 0; j < src.cols(); ++j)
      out(i, j) = scalar_value(src(i, j));
  return out;
}

template<typename Scalar>
inline Eigen::VectorXd to_double_vector(const Vector<Scalar>& src) {
  Eigen::VectorXd out(src.size());
  for (Eigen::Index i = 0; i < src.size(); ++i)
    out(i) = scalar_value(src(i));
  return out;
}

template<typename Scalar>
struct Model {
  Matrix<Scalar> A;
  Matrix<Scalar> Oxx;
  Matrix<Scalar> Oex;
  Matrix<Scalar> Ie;
  Matrix<Scalar> lY;
  Matrix<Scalar> lX;
  Vector<Scalar> tY;
  Vector<Scalar> tX;
  Matrix<Scalar> Gx;
  Matrix<Scalar> Ge;
  Vector<Scalar> a;
  Vector<Scalar> beta0;
  Matrix<Scalar> Psi;
  Matrix<Scalar> d;
  Matrix<Scalar> e;

  unsigned k       = 0;
  unsigned numXis  = 0;
  Eigen::Index deriv_dim = 0;

  Model() = default;

  template<typename OtherScalar>
  explicit Model(const Model<OtherScalar>& other)
    : A(mat_cast<Scalar>(to_double_matrix(other.A), other.deriv_dim)),
      Oxx(mat_cast<Scalar>(to_double_matrix(other.Oxx), other.deriv_dim)),
      Oex(mat_cast<Scalar>(to_double_matrix(other.Oex), other.deriv_dim)),
      Ie(mat_cast<Scalar>(to_double_matrix(other.Ie), other.deriv_dim)),
      lY(mat_cast<Scalar>(to_double_matrix(other.lY), other.deriv_dim)),
      lX(mat_cast<Scalar>(to_double_matrix(other.lX), other.deriv_dim)),
      tY(vec_cast<Scalar>(to_double_vector(other.tY), other.deriv_dim)),
      tX(vec_cast<Scalar>(to_double_vector(other.tX), other.deriv_dim)),
      Gx(mat_cast<Scalar>(to_double_matrix(other.Gx), other.deriv_dim)),
      Ge(mat_cast<Scalar>(to_double_matrix(other.Ge), other.deriv_dim)),
      a(vec_cast<Scalar>(to_double_vector(other.a), other.deriv_dim)),
      beta0(vec_cast<Scalar>(to_double_vector(other.beta0), other.deriv_dim)),
      Psi(mat_cast<Scalar>(to_double_matrix(other.Psi), other.deriv_dim)),
      d(mat_cast<Scalar>(to_double_matrix(other.d), other.deriv_dim)),
      e(mat_cast<Scalar>(to_double_matrix(other.e), other.deriv_dim)),
      k(other.k),
      numXis(other.numXis),
      deriv_dim(other.deriv_dim) {}
};

template<typename Scalar>
inline void set_model_deriv_dim(Model<Scalar>& M, Eigen::Index deriv_dim) {
  M.deriv_dim = deriv_dim;
  if constexpr (ScalarTraits<Scalar>::is_autodiff) {
    auto recast_mat = [&](Matrix<Scalar>& mat) {
      if (mat.size() == 0) return;
      mat = mat_cast<Scalar>(to_double_matrix(mat), deriv_dim);
    };
    auto recast_vec = [&](Vector<Scalar>& vec) {
      if (vec.size() == 0) return;
      vec = vec_cast<Scalar>(to_double_vector(vec), deriv_dim);
    };
    recast_mat(M.A);
    recast_mat(M.Oxx);
    recast_mat(M.Oex);
    recast_mat(M.Ie);
    recast_mat(M.lY);
    recast_mat(M.lX);
    recast_vec(M.tY);
    recast_vec(M.tX);
    recast_mat(M.Gx);
    recast_mat(M.Ge);
    recast_vec(M.a);
    recast_vec(M.beta0);
    recast_mat(M.Psi);
    recast_mat(M.d);
    recast_mat(M.e);
  }
}

struct CompleteData {
  Eigen::MatrixXd V;
  std::vector<std::vector<Eigen::VectorXd>> mean;
  std::vector<std::vector<Eigen::MatrixXd>> cov;
  std::vector<std::vector<double>> tgamma;
  std::vector<Eigen::VectorXi> colidx;
  Eigen::VectorXi n;
  Eigen::VectorXi d;
  int npatterns = 1;
};

struct ObservedData {
  Eigen::MatrixXd V;
  Eigen::VectorXd w;
  std::vector<Eigen::MatrixXd> data;
  std::vector<Eigen::VectorXi> colidx;
  Eigen::VectorXi n;
  int npatterns = 1;
};

template<typename Scalar>
inline Vector<Scalar> make_zvec(unsigned k, unsigned numXis,
                                const Vector<Scalar>& z,
                                Eigen::Index deriv_dim) {
  Vector<Scalar> zVec = zero_vector<Scalar>(numXis, deriv_dim);
  if (k > 0 && z.size() >= static_cast<Eigen::Index>(k)) {
    zVec.head(k) = z.head(k);
  }
  return zVec;
}

template<typename Scalar>
inline Matrix<Scalar> make_oi(unsigned k, unsigned numXis,
                              Eigen::Index deriv_dim) {
  Matrix<Scalar> Oi = zero_matrix<Scalar>(numXis, numXis, deriv_dim);
  for (unsigned i = 0; i < numXis; ++i) {
    double diag_val = (i < k) ? 0.0 : 1.0;
    Oi(i, i) = make_scalar<Scalar>(diag_val, deriv_dim);
  }
  return Oi;
}

template<typename Scalar>
inline Matrix<Scalar> kronecker(const Matrix<Scalar>& A,
                                const Matrix<Scalar>& B) {
  Matrix<Scalar> result(A.rows() * B.rows(), A.cols() * B.cols());
  for (Eigen::Index i = 0; i < A.rows(); ++i) {
    for (Eigen::Index j = 0; j < A.cols(); ++j) {
      result.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) =
        A(i, j) * B;
    }
  }
  return result;
}

template<typename Scalar>
inline Matrix<Scalar> kronecker_vec(const Matrix<Scalar>& A,
                                    const Vector<Scalar>& v) {
  Matrix<Scalar> B(v.size(), 1);
  B.col(0) = v;
  return kronecker(A, B);
}

template<typename Scalar>
inline Vector<Scalar> mu(const Model<Scalar>& M,
                         const Vector<Scalar>& z) {
  using Eigen::numext::log;

  const Vector<Scalar> zVec = make_zvec(M.k, M.numXis, z, M.deriv_dim);
  const Vector<Scalar> beta = M.beta0 + M.A * zVec;

  const Matrix<Scalar> kronZ = kronecker_vec(M.Ie, beta);
  const Matrix<Scalar> Binv =
    (M.Ie - M.Ge - kronZ.transpose() * M.Oex).inverse();

  Vector<Scalar> muX = M.tX + M.lX * beta;
  Vector<Scalar> muY = M.tY +
    M.lY * (Binv * (M.a + M.Gx * beta + kronZ.transpose() * M.Oxx * beta));

  Vector<Scalar> result(muX.size() + muY.size());
  result << muX, muY;
  return result;
}

template<typename Scalar>
inline Matrix<Scalar> sigma(const Model<Scalar>& M,
                            const Vector<Scalar>& z) {
  const Vector<Scalar> zVec = make_zvec(M.k, M.numXis, z, M.deriv_dim);
  const Vector<Scalar> beta = M.beta0 + M.A * zVec;
  const Matrix<Scalar> kronZ = kronecker_vec(M.Ie, beta);
  const Matrix<Scalar> Binv =
    (M.Ie - M.Ge - kronZ.transpose() * M.Oex).inverse();

  Matrix<Scalar> Oi = make_oi<Scalar>(M.k, M.numXis, M.deriv_dim);
  const Matrix<Scalar> Sxx = M.lX * M.A * Oi * M.A.transpose() * M.lX.transpose();
  const Matrix<Scalar> Eta = Binv * (M.Gx * M.A + kronZ.transpose() * M.Oxx * M.A);
  const Matrix<Scalar> Sxy = M.lX * (M.A * Oi * Eta.transpose()) * M.lY.transpose();
  const Matrix<Scalar> Syy = M.lY * Eta * Oi * Eta.transpose() * M.lY.transpose() +
                             M.lY * (Binv * M.Psi * Binv.transpose()) * M.lY.transpose();

  const Eigen::Index p = Sxx.rows();
  const Eigen::Index q = Syy.rows();

  Matrix<Scalar> result(p + q, p + q);
  result.topLeftCorner(p, p) = Sxx;
  result.topRightCorner(p, q) = Sxy;
  result.bottomLeftCorner(q, p) = Sxy.transpose();
  result.bottomRightCorner(q, q) = Syy;

  result += M.d;
  return result;
}

template<typename Scalar>
inline Matrix<Scalar> cholesky_lower(const Matrix<Scalar>& A) {
  const Eigen::Index n = A.rows();
  const Eigen::Index deriv_dim = deriv_dim_of_matrix(A);
  Matrix<Scalar> L = zero_matrix<Scalar>(n, n, deriv_dim);

  for (Eigen::Index i = 0; i < n; ++i) {
    for (Eigen::Index j = 0; j <= i; ++j) {
      Scalar sum = make_scalar<Scalar>(0.0, deriv_dim);
      for (Eigen::Index k = 0; k < j; ++k)
        sum += L(i, k) * L(j, k);

      const Scalar val = A(i, j) - sum;
      if (i == j) {
        L(i, j) = Eigen::numext::sqrt(val);
      } else {
        L(i, j) = val / L(j, j);
      }
    }
  }

  return L;
}

template<typename Scalar>
inline Vector<Scalar>
forward_substitute(const Matrix<Scalar>& L, const Vector<Scalar>& b) {
  const Eigen::Index n = L.rows();
  const Eigen::Index deriv_dim = deriv_dim_of_matrix(L);
  Vector<Scalar> y = zero_vector<Scalar>(n, deriv_dim);

  for (Eigen::Index i = 0; i < n; ++i) {
    Scalar sum = make_scalar<Scalar>(0.0, deriv_dim);
    for (Eigen::Index k = 0; k < i; ++k)
      sum += L(i, k) * y(k);
    y(i) = (b(i) - sum) / L(i, i);
  }
  return y;
}

template<typename Scalar>
inline Vector<Scalar>
backward_substitute(const Matrix<Scalar>& L, const Vector<Scalar>& y) {
  const Eigen::Index n = L.rows();
  const Eigen::Index deriv_dim = deriv_dim_of_matrix(L);
  Vector<Scalar> x = zero_vector<Scalar>(n, deriv_dim);

  for (Eigen::Index i = n; i-- > 0;) {
    Scalar sum = make_scalar<Scalar>(0.0, deriv_dim);
    for (Eigen::Index k = i + 1; k < n; ++k)
      sum += L(k, i) * x(k);
    x(i) = (y(i) - sum) / L(i, i);
  }
  return x;
}

template<typename Scalar>
inline Vector<Scalar>
solve_cholesky(const Matrix<Scalar>& L, const Vector<Scalar>& b) {
  const Vector<Scalar> y = forward_substitute(L, b);
  return backward_substitute(L, y);
}

template<typename Scalar>
inline Matrix<Scalar> cholesky_inverse(const Matrix<Scalar>& L) {
  const Eigen::Index n = L.rows();
  const Eigen::Index deriv_dim = deriv_dim_of_matrix(L);
  Matrix<Scalar> inv = zero_matrix<Scalar>(n, n, deriv_dim);
  for (Eigen::Index i = 0; i < n; ++i) {
    Vector<Scalar> e = zero_vector<Scalar>(n, deriv_dim);
    e(i) = make_scalar<Scalar>(1.0, deriv_dim);
    inv.col(i) = solve_cholesky(L, e);
  }
  return inv;
}

template<typename Scalar>
inline Scalar total_dmvn_weighted(const Vector<Scalar>& mu,
                                  const Matrix<Scalar>& sigma,
                                  const Vector<Scalar>& nu,
                                  const Matrix<Scalar>& S,
                                  double tgamma,
                                  int n,
                                  int d) {
  const double log2pi = std::log(2.0 * M_PI);

  const Eigen::Index deriv_dim = deriv_dim_of_matrix(sigma);
  const Matrix<Scalar> L = cholesky_lower(sigma);

  Scalar log_det = make_scalar<Scalar>(0.0, deriv_dim);
  const Scalar two = make_scalar<Scalar>(2.0, deriv_dim);
  for (Eigen::Index i = 0; i < L.rows(); ++i) {
    log_det += two * Eigen::numext::log(L(i, i));
  }

  const Vector<Scalar> diff = nu - mu;
  const Vector<Scalar> sol = solve_cholesky(L, diff);
  const Scalar tg = make_scalar<Scalar>(tgamma, deriv_dim);
  const Scalar mahal = tg * diff.dot(sol);
  const Matrix<Scalar> sigma_inv = cholesky_inverse(L);
  const Scalar trace_term =
    (sigma_inv.array() * S.array()).sum();

  const Scalar half = make_scalar<Scalar>(0.5, deriv_dim);
  const Scalar log2pi_s = make_scalar<Scalar>(log2pi, deriv_dim);
  const Scalar term = -half * (
    tg * make_scalar<Scalar>(static_cast<double>(d), deriv_dim) * log2pi_s +
    tg * log_det +
    trace_term +
    mahal
  );

  return term;
}

template<typename Scalar>
inline Scalar complete_loglik(const Model<Scalar>& M,
                              const CompleteData& data) {
  Scalar ll = Scalar(0);

  const Matrix<Scalar> V = mat_cast<Scalar>(data.V, M.deriv_dim);

  const auto sum_tgamma = [](const std::vector<double>& tg) {
    double s = 0;
    for (double v : tg) s += v;
    return s;
  };

  for (Eigen::Index j = 0; j < V.rows(); ++j) {
    if (sum_tgamma(data.tgamma[j]) <= std::numeric_limits<double>::min()) {
      continue;
    }

    const Vector<Scalar> z = V.row(j).transpose();
    const Vector<Scalar> mu_vec = mu(M, z);
    const Matrix<Scalar> Sig = sigma(M, z);

    for (int i = 0; i < data.npatterns; ++i) {
      const double tg = data.tgamma[j][i];
      if (tg <= std::numeric_limits<double>::min()) continue;

      const Eigen::VectorXi& cols = data.colidx[i];
      Vector<Scalar> mu_sub = zero_vector<Scalar>(cols.size(), M.deriv_dim);
      Matrix<Scalar> Sig_sub = zero_matrix<Scalar>(cols.size(), cols.size(), M.deriv_dim);

      for (Eigen::Index r = 0; r < cols.size(); ++r) {
        mu_sub(r) = mu_vec(cols[r]);
        for (Eigen::Index c = 0; c < cols.size(); ++c) {
          Sig_sub(r, c) = Sig(cols[r], cols[c]);
        }
      }

      const Vector<Scalar> nu = vec_cast<Scalar>(data.mean[j][i], M.deriv_dim);
      const Matrix<Scalar> S  = mat_cast<Scalar>(data.cov[j][i], M.deriv_dim);
      ll += total_dmvn_weighted(
        mu_sub, Sig_sub, nu, S, tg, data.n[i], data.d[i]
      );
    }
  }

  return ll;
}

template<typename Scalar>
inline Scalar observed_loglik(const Model<Scalar>& M,
                              const ObservedData& data) {
  const Matrix<Scalar> V = mat_cast<Scalar>(data.V, M.deriv_dim);
  const Vector<Scalar> w = vec_cast<Scalar>(data.w, M.deriv_dim);
  const int total_n = data.n.sum();

  Vector<Scalar> density = zero_vector<Scalar>(total_n, M.deriv_dim);

  for (Eigen::Index i = 0; i < V.rows(); ++i) {
    if (w(i) <= make_scalar<Scalar>(std::numeric_limits<double>::min(), M.deriv_dim))
      continue;

    const Vector<Scalar> z = V.row(i).transpose();
    const Vector<Scalar> mu_vec = mu(M, z);
    const Matrix<Scalar> Sig = sigma(M, z);

    int idx_offset = 0;
    for (int j = 0; j < data.npatterns; ++j) {
      const Eigen::VectorXi& cols = data.colidx[j];
      Vector<Scalar> mu_sub = zero_vector<Scalar>(cols.size(), M.deriv_dim);
      Matrix<Scalar> Sig_sub = zero_matrix<Scalar>(cols.size(), cols.size(), M.deriv_dim);

      for (Eigen::Index r = 0; r < cols.size(); ++r) {
        mu_sub(r) = mu_vec(cols[r]);
        for (Eigen::Index c = 0; c < cols.size(); ++c) {
          Sig_sub(r, c) = Sig(cols[r], cols[c]);
        }
      }

      const Eigen::MatrixXd& X = data.data[j];
      const Matrix<Scalar> L = cholesky_lower(Sig_sub);

      Scalar log_det = make_scalar<Scalar>(0.0, M.deriv_dim);
      const Scalar two = make_scalar<Scalar>(2.0, M.deriv_dim);
      for (Eigen::Index r = 0; r < L.rows(); ++r) {
        log_det += two * Eigen::numext::log(L(r, r));
      }

      const int d = cols.size();
      const Scalar half = make_scalar<Scalar>(0.5, M.deriv_dim);
      const Scalar norm_const =
        -half * make_scalar<Scalar>(static_cast<double>(d), M.deriv_dim) *
        Eigen::numext::log(make_scalar<Scalar>(2.0 * M_PI, M.deriv_dim));
      const Scalar det_const = -half * log_det;

      for (Eigen::Index row = 0; row < X.rows(); ++row) {
        Vector<Scalar> xvec = zero_vector<Scalar>(cols.size(), M.deriv_dim);
        for (Eigen::Index c = 0; c < cols.size(); ++c) {
          xvec(c) = make_scalar<Scalar>(X(row, c), M.deriv_dim);
        }
        Vector<Scalar> diff = xvec - mu_sub;
        const Vector<Scalar> sol = solve_cholesky(L, diff);
        const Scalar quad = -half * diff.dot(sol);
        const Scalar log_pdf = norm_const + det_const + quad;
        density(idx_offset + row) +=
          Eigen::numext::exp(log_pdf) * w(i);
      }

      idx_offset += X.rows();
    }
  }

  Scalar loglik = make_scalar<Scalar>(0.0, M.deriv_dim);
  for (Eigen::Index r = 0; r < density.size(); ++r) {
    loglik += Eigen::numext::log(density(r));
  }
  return loglik;
}

} // namespace lms
} // namespace modsem

namespace Eigen {

template<>
struct NumTraits<modsem::lms::Dual2> : NumTraits<double> {
  typedef modsem::lms::Dual2 Real;
  typedef modsem::lms::Dual2 NonInteger;
  typedef modsem::lms::Dual2 Nested;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 3,
    AddCost = 3,
    MulCost = 3
  };
};

namespace numext {
inline modsem::lms::Dual2 sqrt(const modsem::lms::Dual2& x) { return modsem::lms::sqrt(x); }
inline modsem::lms::Dual2 log(const modsem::lms::Dual2& x) { return modsem::lms::log(x); }
inline modsem::lms::Dual2 exp(const modsem::lms::Dual2& x) { return modsem::lms::exp(x); }
inline modsem::lms::Dual2 abs(const modsem::lms::Dual2& x) { return modsem::lms::abs(x); }
inline bool isfinite(const modsem::lms::Dual2& x) { return std::isfinite(x.val); }
inline bool isnan(const modsem::lms::Dual2& x) { return std::isnan(x.val); }
inline bool isinf(const modsem::lms::Dual2& x) { return std::isinf(x.val); }
}

} // namespace Eigen

#endif
