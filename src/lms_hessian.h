#ifndef MODSEM_LMS_HESSIAN_H
#define MODSEM_LMS_HESSIAN_H

// Matrix-level forward-over-reverse differentiation for the LMS Hessian.
// Values and reverse adjoints are cached once per quadrature node. Each column
// then differentiates the reverse pass analytically; no parameter perturbations
// or repeated covariance factorizations are needed. Memory is O(graph size + p^2),
// rather than storing a p by p Hessian for every matrix element.
namespace lms_hessian {

struct Dual {
  arma::mat v, d;
  Dual(const arma::mat& value, const arma::mat& tangent) : v(value), d(tangent) {}
  Dual t() const { return Dual(v.t(), d.t()); }
};
inline Dual operator+(const Dual& a, const Dual& b) {
  return Dual(a.v + b.v, a.d + b.d);
}
inline Dual operator-(const Dual& a, const Dual& b) {
  return Dual(a.v - b.v, a.d - b.d);
}
inline Dual operator-(const Dual& a) { return Dual(-a.v, -a.d); }
inline Dual operator*(const Dual& a, const Dual& b) {
  return Dual(a.v * b.v, a.d * b.v + a.v * b.d);
}

enum class Op { Input, Add, Subtract, Multiply, Transpose, Inverse, Pinverse,
                Kron, JoinRows, JoinCols };

struct Node {
  Op op;
  int a, b;
  bool variable = false, active = false, barActive = false;
  arma::mat v, dot, bar, barDot;
  Node(Op op_, int a_, int b_, const arma::mat& value) :
    op(op_), a(a_), b(b_), v(value),
    dot(value.n_rows, value.n_cols, arma::fill::zeros),
    bar(value.n_rows, value.n_cols, arma::fill::zeros),
    barDot(value.n_rows, value.n_cols, arma::fill::zeros) {}
};

class Tape {
public:
  std::vector<Node> nodes;

  int input(const arma::mat& v) { return append(Op::Input, -1, -1, v); }
  int add(int a, int b) { return append(Op::Add, a, b, nodes[a].v + nodes[b].v); }
  int sub(int a, int b) { return append(Op::Subtract, a, b, nodes[a].v - nodes[b].v); }
  int mul(int a, int b) { return append(Op::Multiply, a, b, nodes[a].v * nodes[b].v); }
  int trans(int a) { return append(Op::Transpose, a, -1, nodes[a].v.t()); }
  int inv(int a) { return append(Op::Inverse, a, -1, arma::inv(nodes[a].v)); }
  int pinv(int a) { return append(Op::Pinverse, a, -1, arma::pinv(nodes[a].v)); }
  // The left Kronecker factor is the constant Ieta.
  int kron(int a, int b) { return append(Op::Kron, a, b, arma::kron(nodes[a].v, nodes[b].v)); }
  int rows(int a, int b) { return append(Op::JoinRows, a, b, arma::join_rows(nodes[a].v, nodes[b].v)); }
  int cols(int a, int b) { return append(Op::JoinCols, a, b, arma::join_cols(nodes[a].v, nodes[b].v)); }

  void forward() {
    for (auto& y : nodes) {
      if (y.op == Op::Input) continue;
      const Node& a = nodes[y.a];
      y.active = a.active || (y.b >= 0 && nodes[y.b].active);
      if (!y.active) continue;
      switch (y.op) {
      case Op::Add: y.dot = a.dot + nodes[y.b].dot; break;
      case Op::Subtract: y.dot = a.dot - nodes[y.b].dot; break;
      case Op::Multiply: {
        const Node& b = nodes[y.b];
        if (a.active) y.dot = a.dot * b.v;
        if (b.active) y.dot += a.v * b.dot;
        break;
      }
      case Op::Transpose: y.dot = a.dot.t(); break;
      case Op::Inverse: y.dot = -y.v * a.dot * y.v; break;
      case Op::Pinverse: {
        const arma::mat I = arma::eye<arma::mat>(a.v.n_rows, a.v.n_cols);
        y.dot = -y.v * a.dot * y.v
          + y.v * y.v.t() * a.dot.t() * (I - a.v * y.v)
          + (I - y.v * a.v) * a.dot.t() * y.v.t() * y.v;
        break;
      }
      case Op::Kron: y.dot = arma::kron(a.v, nodes[y.b].dot); break;
      case Op::JoinRows: y.dot = arma::join_rows(a.dot, nodes[y.b].dot); break;
      case Op::JoinCols: y.dot = arma::join_cols(a.dot, nodes[y.b].dot); break;
      default: break;
      }
    }
  }

  // First reverse pass computes bars; subsequent passes compute their tangents.
  void reverse(bool tangent) {
    for (std::size_t j = nodes.size(); j-- > 0;) {
      Node& y = nodes[j];
      if (y.op == Op::Input || !y.variable) continue;
      Node& a = nodes[y.a];
      if (tangent && !y.barActive &&
          !(y.op == Op::Multiply && (a.active || nodes[y.b].active)) &&
          !((y.op == Op::Inverse || y.op == Op::Pinverse) && y.active)) continue;
      arma::mat& ab = tangent ? a.barDot : a.bar;
      const arma::mat& yb = tangent ? y.barDot : y.bar;
      switch (y.op) {
      case Op::Add:
        if (a.variable) ab += yb;
        if (nodes[y.b].variable) (tangent ? nodes[y.b].barDot : nodes[y.b].bar) += yb;
        break;
      case Op::Subtract:
        if (a.variable) ab += yb;
        if (nodes[y.b].variable) (tangent ? nodes[y.b].barDot : nodes[y.b].bar) -= yb;
        break;
      case Op::Multiply: {
        Node& b = nodes[y.b];
        if (tangent) {
          if (y.barActive) {
            if (a.variable) ab += y.barDot * b.v.t();
            if (b.variable) b.barDot += a.v.t() * y.barDot;
          }
          if (a.variable && b.active) ab += y.bar * b.dot.t();
          if (b.variable && a.active) b.barDot += a.dot.t() * y.bar;
        } else {
          if (a.variable) ab += y.bar * b.v.t();
          if (b.variable) b.bar += a.v.t() * y.bar;
        }
        break;
      }
      case Op::Transpose: ab += yb.t(); break;
      case Op::Inverse:
        if (tangent) {
          if (y.active) ab -= y.dot.t() * y.bar * y.v.t()
            + y.v.t() * y.bar * y.dot.t();
          if (y.barActive) ab -= y.v.t() * y.barDot * y.v.t();
        } else ab -= y.v.t() * y.bar * y.v.t();
        break;
      case Op::Pinverse: {
        // Moore-Penrose differential on a constant-rank neighbourhood. The
        // null-space terms matter for composite models with singular W'TW.
        const arma::mat I = arma::eye<arma::mat>(a.v.n_rows, a.v.n_cols);
        if (tangent) {
          const Dual A(a.v, a.dot), R(y.v, y.dot), G(y.bar, y.barDot);
          const Dual Id(I, arma::zeros<arma::mat>(I.n_rows, I.n_cols));
          const Dual result = -R.t() * G * R.t()
            + (R * R.t()) * G * (Id - A * R).t()
            + (Id - R * A).t() * G * (R.t() * R);
          ab += result.d;
        } else {
          ab += -y.v.t() * y.bar * y.v.t()
            + (y.v * y.v.t()) * y.bar * (I - a.v * y.v).t()
            + (I - y.v * a.v).t() * y.bar * (y.v.t() * y.v);
        }
        break;
      }
      case Op::Kron: {
        Node& b = nodes[y.b];
        arma::mat& bb = tangent ? b.barDot : b.bar;
        for (arma::uword c = 0; c < a.v.n_cols; ++c)
          for (arma::uword r = 0; r < a.v.n_rows; ++r)
            if (a.v(r, c) != 0.0)
              bb += a.v(r, c) * yb.submat(r * b.v.n_rows, c * b.v.n_cols,
                (r + 1) * b.v.n_rows - 1, (c + 1) * b.v.n_cols - 1);
        break;
      }
      case Op::JoinRows: {
        Node& b = nodes[y.b];
        if (a.variable && a.v.n_cols) ab += yb.head_cols(a.v.n_cols);
        if (b.variable && b.v.n_cols) (tangent ? b.barDot : b.bar) += yb.tail_cols(b.v.n_cols);
        break;
      }
      case Op::JoinCols: {
        Node& b = nodes[y.b];
        if (a.variable && a.v.n_rows) ab += yb.head_rows(a.v.n_rows);
        if (b.variable && b.v.n_rows) (tangent ? b.barDot : b.bar) += yb.tail_rows(b.v.n_rows);
        break;
      }
      default: break;
      }
      if (tangent) {
        a.barActive = a.barActive || y.barActive ||
          (y.op == Op::Multiply && nodes[y.b].active) ||
          ((y.op == Op::Inverse || y.op == Op::Pinverse) && y.active);
        if (y.b >= 0) nodes[y.b].barActive = nodes[y.b].barActive ||
          y.barActive || (y.op == Op::Multiply && a.active);
      }
    }
  }

  void resetDirection() {
    for (auto& node : nodes) {
      node.dot.zeros();
      node.barDot.zeros();
      node.active = node.barActive = false;
    }
  }

private:
  int append(Op op, int a, int b, const arma::mat& v) {
    nodes.emplace_back(op, a, b, v);
    return static_cast<int>(nodes.size() - 1);
  }
};

struct Moments {
  Tape tape;
  std::vector<int> blocks;
  int mu, sigma;

  Moments(const LMSModel& M, const arma::vec& z,
          const arma::uvec& requestedBlocks) : blocks(18, -1) {
    Tape& t = tape;
    t.nodes.reserve(100);
    const int A = blocks[6] = t.input(M.A);
    const int Oxx = blocks[12] = t.input(M.Oxx);
    const int Oex = blocks[13] = t.input(M.Oex);
    const int lX = blocks[0] = t.input(M.lX);
    const int tau = blocks[2] = t.input(M.tX);
    const int Gx = blocks[10] = t.input(M.Gx);
    const int Ge = blocks[11] = t.input(M.Ge);
    const int a = blocks[8] = t.input(M.a);
    const int beta = blocks[9] = t.input(M.beta0);
    const int psi = blocks[7] = t.input(M.Psi);
    const int cov = blocks[17] = t.input(M.covZetaXi);
    const int d = blocks[4] = t.input(M.d);
    const int Ie = t.input(M.Ie), Oi = t.input(M.Oi);
    const int zv = t.input(makeZvec(M.k, M.numXis, z));
    const int proj = t.mul(cov, t.trans(t.inv(A)));
    const int orth = t.sub(psi, t.mul(proj, t.trans(proj)));
    const int u = t.add(beta, t.mul(A, zv));
    const int K = t.kron(Ie, u);
    const int KO = t.mul(t.trans(K), Oxx);
    const int B = t.inv(t.sub(t.sub(Ie, Ge), t.mul(t.trans(K), Oex)));
    const int r = t.add(t.add(a, t.mul(Gx, u)),
                        t.add(t.mul(KO, u), t.mul(proj, zv)));
    const int v = t.mul(B, r);
    const int C = t.add(t.mul(t.add(Gx, KO), A), proj);
    const int Eta = t.mul(B, C);
    const int AO = t.mul(A, Oi);
    const int varXi = t.mul(AO, t.trans(A));
    const int cross = t.mul(AO, t.trans(Eta));
    const int varEta = t.add(t.mul(t.mul(Eta, Oi), t.trans(Eta)),
                             t.mul(t.mul(B, orth), t.trans(B)));
    const int V = t.cols(t.rows(varXi, cross), t.rows(t.trans(cross), varEta));
    const int xieta = t.cols(u, v);
    int L = lX, latentCov = V, residual = d;
    if (M.hasComposites) {
      const int W = blocks[14] = t.input(M.W);
      const int T = blocks[15] = t.input(M.T);
      const int TW = t.mul(T, W);
      const int H = t.mul(t.trans(W), TW);
      L = t.add(lX, t.mul(TW, t.pinv(H)));
      latentCov = t.sub(V, H);
      residual = t.add(d, T);
    }
    mu = t.add(tau, t.mul(L, xieta));
    sigma = t.add(t.mul(t.mul(L, latentCov), t.trans(L)), residual);
    // Constants and model matrices without requested parameters never need
    // reverse adjoints, but their values still enter other derivatives.
    for (unsigned block : requestedBlocks) {
      const int id = blocks.at(block);
      if (id >= 0) t.nodes[id].variable = true;
    }
    for (auto& node : t.nodes)
      if (node.op != Op::Input)
        node.variable = t.nodes[node.a].variable ||
          (node.b >= 0 && t.nodes[node.b].variable);
  }

  double derivative(unsigned block, unsigned row, unsigned col,
                    bool symmetric, bool tangent) const {
    const int id = blocks.at(block);
    if (id < 0) return 0.0;
    const Node& node = tape.nodes[id];
    const arma::mat& bar = tangent ? node.barDot : node.bar;
    return bar(row, col) + ((symmetric && row != col) ? bar(col, row) : 0.0);
  }
};

// Gaussian score and its directional derivative with the posterior sufficient
// statistics held fixed. S is the weighted scatter, not a normalized covariance.
struct Pattern {
  arma::uvec idx;
  arma::mat inverse, sigmaBar;
  arma::vec inverseDiff;
  double weight;

  bool initialize(Moments& m, const arma::uvec& indices, const arma::vec& mean,
                  const arma::mat& cov, double tg, double& value) {
    idx = indices;
    weight = tg;
    Node& mu = m.tape.nodes[m.mu];
    Node& sigma = m.tape.nodes[m.sigma];
    const arma::mat sig = sigma.v.submat(idx, idx);
    arma::mat L;
    if (!sig.is_finite() || !arma::chol(L, sig, "lower")) return false;
    inverse = arma::solve(arma::trimatu(L.t()),
      arma::solve(arma::trimatl(L), arma::eye<arma::mat>(idx.n_elem, idx.n_elem),
                  arma::solve_opts::fast), arma::solve_opts::fast);
    const arma::vec diff = mean - mu.v.elem(idx);
    inverseDiff = inverse * diff;
    const arma::mat scatter = cov + tg * diff * diff.t();
    sigmaBar = 0.5 * (inverse * scatter * inverse - tg * inverse);
    mu.bar.elem(idx) += tg * inverseDiff;
    sigma.bar.submat(idx, idx) += sigmaBar;
    value -= 0.5 * (tg * (idx.n_elem * std::log(2.0 * M_PI)
      + 2.0 * arma::sum(arma::log(L.diag()))) + arma::accu(inverse % scatter.t()));
    return true;
  }

  void seedDirection(Moments& m) const {
    Node& mu = m.tape.nodes[m.mu];
    Node& sigma = m.tape.nodes[m.sigma];
    const arma::vec dm = mu.dot.elem(idx);
    const arma::mat ds = sigma.dot.submat(idx, idx);
    if (!mu.active && !sigma.active) return;
    mu.barActive = sigma.barActive = true;
    const arma::vec inverseDm = inverse * dm;
    mu.barDot.elem(idx) -= weight * inverseDm;
    sigma.barDot.submat(idx, idx) -= 0.5 * weight * (
      inverseDm * inverseDiff.t() + inverseDiff * inverseDm.t());
    if (sigma.active) {
      // If Q = Sigma^-1 and G = dL/dSigma, then
      // dG = -Q dSigma G - G dSigma Q - (weight/2) Q dSigma Q
      //      - (weight/2) [(Q dmu)(Q diff)' + (Q diff)(Q dmu)'].
      const arma::mat inverseDs = inverse * ds;
      mu.barDot.elem(idx) -= weight * inverseDs * inverseDiff;
      sigma.barDot.submat(idx, idx) -= inverseDs * sigmaBar
        + sigmaBar * inverseDs.t() + 0.5 * weight * inverseDs * inverse;
    }
  }
};

} // namespace lms_hessian
#endif
