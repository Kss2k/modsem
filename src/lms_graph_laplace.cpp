#include <RcppArmadillo.h>
#include <stdexcept>
// [[Rcpp::depends(RcppArmadillo)]]
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Jet4 {
  int d;
  double v;
  arma::vec g;
  arma::mat h;
  std::vector<double> t;
  std::vector<double> q;

  Jet4(const int dimension = 0, const double value = 0.0) :
    d(dimension), v(value), g(dimension, arma::fill::zeros),
    h(dimension, dimension, arma::fill::zeros),
    t(static_cast<std::size_t>(dimension) * dimension * dimension, 0.0),
    q(static_cast<std::size_t>(dimension) * dimension * dimension * dimension,
      0.0) {}

  inline std::size_t i3(int i, int j, int k) const {
    return (static_cast<std::size_t>(i) * d + j) * d + k;
  }
  inline std::size_t i4(int i, int j, int k, int l) const {
    return ((static_cast<std::size_t>(i) * d + j) * d + k) * d + l;
  }
};

Jet4 constantJet(const int d, const double value) { return Jet4(d, value); }

Jet4 variableJet(const int d, const int index, const double value) {
  Jet4 out(d, value);
  out.g(index) = 1.0;
  return out;
}

Jet4 addJet(const Jet4& a, const Jet4& b) {
  if(a.d!=b.d) throw std::runtime_error("Laplace jet dimension mismatch in addition.");
  if(a.t.size()!=b.t.size()||a.q.size()!=b.q.size())
    throw std::runtime_error("Laplace jet tensor mismatch in addition.");
  Jet4 out(a.d, a.v + b.v);
  out.g = a.g + b.g;
  out.h = a.h + b.h;
  for (std::size_t i = 0; i < out.t.size(); ++i) out.t.at(i) = a.t.at(i) + b.t.at(i);
  for (std::size_t i = 0; i < out.q.size(); ++i) out.q.at(i) = a.q.at(i) + b.q.at(i);
  return out;
}

Jet4 scaleJet(const Jet4& a, const double scale) {
  Jet4 out(a.d, scale * a.v);
  out.g = scale * a.g;
  out.h = scale * a.h;
  if(out.t.size()!=a.t.size()||out.q.size()!=a.q.size())
    throw std::runtime_error("Laplace jet tensor mismatch in scaling.");
  for (std::size_t i = 0; i < out.t.size(); ++i) out.t.at(i) = scale * a.t.at(i);
  for (std::size_t i = 0; i < out.q.size(); ++i) out.q.at(i) = scale * a.q.at(i);
  return out;
}

Jet4 multiplyJet(const Jet4& a, const Jet4& b) {
  if(a.d!=b.d) throw std::runtime_error("Laplace jet dimension mismatch in multiplication.");
  const int d = a.d;
  Jet4 out(d, a.v * b.v);
  for (int i = 0; i < d; ++i) {
    out.g(i) = a.g(i) * b.v + a.v * b.g(i);
    for (int j = 0; j < d; ++j) {
      out.h(i, j) = a.h(i, j) * b.v + a.g(i) * b.g(j) +
        a.g(j) * b.g(i) + a.v * b.h(i, j);
      for (int k = 0; k < d; ++k) {
        const std::size_t ijk = out.i3(i, j, k);
        out.t.at(ijk) = a.t.at(ijk) * b.v +
          a.h(i, j) * b.g(k) + a.h(i, k) * b.g(j) +
          a.h(j, k) * b.g(i) +
          a.g(i) * b.h(j, k) + a.g(j) * b.h(i, k) +
          a.g(k) * b.h(i, j) + a.v * b.t.at(ijk);
        for (int l = 0; l < d; ++l) {
          const std::size_t ijkl = out.i4(i, j, k, l);
          out.q.at(ijkl) = a.q.at(ijkl) * b.v + a.v * b.q.at(ijkl) +
            a.t.at(out.i3(i,j,k)) * b.g(l) + a.t.at(out.i3(i,j,l)) * b.g(k) +
            a.t.at(out.i3(i,k,l)) * b.g(j) + a.t.at(out.i3(j,k,l)) * b.g(i) +
            a.g(i) * b.t.at(out.i3(j,k,l)) + a.g(j) * b.t.at(out.i3(i,k,l)) +
            a.g(k) * b.t.at(out.i3(i,j,l)) + a.g(l) * b.t.at(out.i3(i,j,k)) +
            a.h(i,j) * b.h(k,l) + a.h(i,k) * b.h(j,l) +
            a.h(i,l) * b.h(j,k) + a.h(j,k) * b.h(i,l) +
            a.h(j,l) * b.h(i,k) + a.h(k,l) * b.h(i,j);
        }
      }
    }
  }
  return out;
}

Jet4 composeJet(const Jet4& u, const double value, const double d1,
                const double d2, const double d3, const double d4) {
  const int d = u.d;
  Jet4 out(d, value);
  for (int i = 0; i < d; ++i) {
    out.g(i) = d1 * u.g(i);
    for (int j = 0; j < d; ++j) {
      out.h(i,j) = d2 * u.g(i) * u.g(j) + d1 * u.h(i,j);
      for (int k = 0; k < d; ++k) {
        const std::size_t ijk = out.i3(i,j,k);
        out.t.at(ijk) = d3 * u.g(i) * u.g(j) * u.g(k) +
          d2 * (u.h(i,j) * u.g(k) + u.h(i,k) * u.g(j) +
                u.h(j,k) * u.g(i)) + d1 * u.t.at(ijk);
        for (int l = 0; l < d; ++l) {
          const std::size_t ijkl = out.i4(i,j,k,l);
          const double hgg = u.h(i,j)*u.g(k)*u.g(l) +
            u.h(i,k)*u.g(j)*u.g(l) + u.h(i,l)*u.g(j)*u.g(k) +
            u.h(j,k)*u.g(i)*u.g(l) + u.h(j,l)*u.g(i)*u.g(k) +
            u.h(k,l)*u.g(i)*u.g(j);
          const double hh = u.h(i,j)*u.h(k,l) + u.h(i,k)*u.h(j,l) +
            u.h(i,l)*u.h(j,k);
          const double tg = u.t.at(u.i3(i,j,k))*u.g(l) +
            u.t.at(u.i3(i,j,l))*u.g(k) + u.t.at(u.i3(i,k,l))*u.g(j) +
            u.t.at(u.i3(j,k,l))*u.g(i);
          out.q.at(ijkl) = d4*u.g(i)*u.g(j)*u.g(k)*u.g(l) +
            d3*hgg + d2*(hh + tg) + d1*u.q.at(ijkl);
        }
      }
    }
  }
  return out;
}

Jet4 logJet(const Jet4& u) {
  const double inverse = 1.0 / u.v;
  return composeJet(u, std::log(u.v), inverse, -inverse*inverse,
                    2.0*std::pow(inverse,3), -6.0*std::pow(inverse,4));
}

Jet4 logisticCdfJet(const Jet4& u) {
  const double s = 1.0 / (1.0 + std::exp(-u.v));
  const double d1 = s * (1.0 - s);
  const double d2 = d1 * (1.0 - 2.0*s);
  const double d3 = d1 * (1.0 - 6.0*s + 6.0*s*s);
  const double d4 = d1 * (1.0 - 14.0*s + 36.0*s*s - 24.0*s*s*s);
  return composeJet(u, s, d1, d2, d3, d4);
}

Jet4 probitCdfJet(const Jet4& u) {
  const double pdf = R::dnorm(u.v, 0.0, 1.0, false);
  return composeJet(u, R::pnorm(u.v, 0.0, 1.0, true, false), pdf,
    -u.v*pdf, (u.v*u.v - 1.0)*pdf,
    (-u.v*u.v*u.v + 3.0*u.v)*pdf);
}

Jet4 cdfJet(const Jet4& u, const bool logistic) {
  return logistic ? logisticCdfJet(u) : probitCdfJet(u);
}

struct Term { int first; int second; double coefficient; };
struct Equation {
  std::vector<Term> xiLinear, etaLinear, xiXi, xiEta;
};

std::vector<Equation> graphPlan(const arma::mat& gx, const arma::mat& ge,
                                const arma::mat& oxx, const arma::mat& oex,
                                const int nx, const int ne) {
  std::vector<Equation> plan(ne);
  for (int equation = 0; equation < ne; ++equation) {
    const int firstRow = equation * nx;
    for (int x = 0; x < nx; ++x) {
      if (gx(equation,x) != 0.0) plan[equation].xiLinear.push_back({x,-1,gx(equation,x)});
      for (int z = 0; z < nx; ++z)
        if (oxx(firstRow+x,z) != 0.0)
          plan[equation].xiXi.push_back({x,z,oxx(firstRow+x,z)});
      for (int e = 0; e < ne; ++e)
        if (oex(firstRow+x,e) != 0.0)
          plan[equation].xiEta.push_back({x,e,oex(firstRow+x,e)});
    }
    for (int e = 0; e < ne; ++e)
      if (ge(equation,e) != 0.0)
        plan[equation].etaLinear.push_back({e,-1,ge(equation,e)});
  }
  return plan;
}

struct Observation { arma::rowvec values; std::vector<int> columns; };

struct LaplaceModel {
  int nx, ne, d;
  arma::mat chol, gx, ge, oxx, oex, lambda, theta;
  arma::vec beta0, alpha, tau;
  std::vector<Equation> plan;
  std::vector<bool> ordered;
  std::vector<arma::rowvec> thresholds;
  bool logistic;

  LaplaceModel(const Rcpp::List& matrices, const int numXis,
               const int numEtas, const Rcpp::IntegerVector& orderedIndex,
               const bool useLogistic) : nx(numXis), ne(numEtas),
    d(numXis + numEtas),
    gx(Rcpp::as<arma::mat>(matrices["gammaXi"])),
    ge(Rcpp::as<arma::mat>(matrices["gammaEta"])),
    oxx(Rcpp::as<arma::mat>(matrices["omegaXiXi"])),
    oex(Rcpp::as<arma::mat>(matrices["omegaEtaXi"])),
    lambda(Rcpp::as<arma::mat>(matrices["lambdaX"])),
    theta(Rcpp::as<arma::mat>(matrices["thetaDelta"])),
    beta0(Rcpp::as<arma::vec>(matrices["beta0"])),
    alpha(Rcpp::as<arma::vec>(matrices["alpha"])),
    tau(Rcpp::as<arma::vec>(matrices["tauX"])), logistic(useLogistic) {
    const arma::mat A=Rcpp::as<arma::mat>(matrices["A"]);
    const arma::mat cross=Rcpp::as<arma::mat>(matrices["covZetaXi"]);
    const arma::mat psi=Rcpp::as<arma::mat>(matrices["psi"]);
    arma::mat covariance(d,d,arma::fill::zeros);
    if(nx) covariance.submat(0,0,nx-1,nx-1)=A*A.t();
    if(ne) covariance.submat(nx,nx,d-1,d-1)=psi;
    if(nx&&ne) {
      covariance.submat(nx,0,d-1,nx-1)=cross;
      covariance.submat(0,nx,nx-1,d-1)=cross.t();
    }
    covariance=0.5*(covariance+covariance.t());
    if(!arma::chol(chol,covariance))
      Rcpp::stop("The unified latent covariance matrix is not positive definite.");
    plan=graphPlan(gx,ge,oxx,oex,nx,ne);
    ordered.assign(lambda.n_rows,false);
    for(int index:orderedIndex)
      if(index>=0&&static_cast<arma::uword>(index)<lambda.n_rows)
        ordered[index]=true;
    thresholds.resize(lambda.n_rows);
    const arma::mat thresholdMatrix=Rcpp::as<arma::mat>(matrices["thresholds"]);
    for(arma::uword column=0;column<lambda.n_rows;++column) if(ordered[column]) {
      const arma::uvec use=arma::find_finite(thresholdMatrix.row(column).t());
      thresholds[column].set_size(use.n_elem);
      for(arma::uword k=0;k<use.n_elem;++k)
        thresholds[column](k)=thresholdMatrix(column,use(k));
    }
  }
};

std::vector<Observation> observationsFromR(const Rcpp::List& dataR,
                                           const Rcpp::List& colidxR) {
  std::vector<Observation> out;
  for (int p = 0; p < dataR.size(); ++p) {
    const arma::mat values = Rcpp::as<arma::mat>(dataR[p]);
    const Rcpp::IntegerVector columnsR = colidxR[p];
    const std::vector<int> columns(columnsR.begin(), columnsR.end());
    for (arma::uword i = 0; i < values.n_rows; ++i)
      out.push_back({values.row(i), columns});
  }
  return out;
}

struct LaplaceResult {
  double first;
  double second;
  double correction;
  bool adjusted;
  arma::vec gradient;
  arma::mat hessian;
};

LaplaceResult observationLaplace(
    const LaplaceModel& model, const arma::rowvec& mode,
    const Observation& observation, const double curvatureFloor,
    const double correctionFloor) {
  const int nx=model.nx, ne=model.ne, d=model.d;
  std::vector<Jet4> z, innovations(d), states(d);
  for(int j=0;j<d;++j) states[j]=constantJet(d,0.0);
  for (int j=0;j<d;++j) z.push_back(variableJet(d,j,mode(j)));
  for (int j=0;j<d;++j) {
    innovations[j]=constantJet(d,0.0);
    for (int k=0;k<d;++k)
      innovations[j]=addJet(innovations[j],scaleJet(z[k],model.chol(k,j)));
  }
  for (int j=0;j<nx;++j) states[j]=addJet(innovations[j],constantJet(d,model.beta0(j)));
  for (int j=0;j<ne;++j) {
    Jet4 value=addJet(innovations[nx+j],constantJet(d,model.alpha(j)));
    for(const Term& x:model.plan[j].xiLinear) value=addJet(value,scaleJet(states[x.first],x.coefficient));
    for(const Term& e:model.plan[j].etaLinear) value=addJet(value,scaleJet(states[nx+e.first],e.coefficient));
    for(const Term& x:model.plan[j].xiXi) value=addJet(value,scaleJet(multiplyJet(states[x.first],states[x.second]),x.coefficient));
    for(const Term& x:model.plan[j].xiEta) value=addJet(value,scaleJet(multiplyJet(states[x.first],states[nx+x.second]),x.coefficient));
    states[nx+j]=value;
  }
  Jet4 f=constantJet(d,0.5*d*std::log(2.0*M_PI));
  for(int j=0;j<d;++j) f=addJet(f,scaleJet(multiplyJet(z[j],z[j]),0.5));
  for(arma::uword jj=0;jj<observation.values.n_elem;++jj) {
    const int column=observation.columns[jj];
    Jet4 mu=constantJet(d,model.tau(column));
    for(int latent=0;latent<d;++latent)
      if(model.lambda(column,latent)!=0.0)
        mu=addJet(mu,scaleJet(states[latent],model.lambda(column,latent)));
    if(model.ordered[column]) {
      const int code=static_cast<int>(observation.values(jj));
      const int categories=model.thresholds[column].n_elem+1;
      Jet4 probability=constantJet(d,1.0);
      if(code<categories) probability=cdfJet(addJet(constantJet(d,model.thresholds[column](code-1)),scaleJet(mu,-1.0)),model.logistic);
      if(code>1) probability=addJet(probability,scaleJet(cdfJet(addJet(constantJet(d,model.thresholds[column](code-2)),scaleJet(mu,-1.0)),model.logistic),-1.0));
      if(!(probability.v>0.0) || !std::isfinite(probability.v))
        return {NA_REAL,NA_REAL,NA_REAL,true,
          arma::vec(d,arma::fill::value(NA_REAL)),
          arma::mat(d,d,arma::fill::value(NA_REAL))};
      f=addJet(f,scaleJet(logJet(probability),-1.0));
    } else {
      const double variance=model.theta(column,column);
      Jet4 residual=addJet(constantJet(d,observation.values(jj)),scaleJet(mu,-1.0));
      f=addJet(f,constantJet(d,0.5*std::log(2.0*M_PI*variance)));
      f=addJet(f,scaleJet(multiplyJet(residual,residual),0.5/variance));
    }
  }
  arma::vec eigval; arma::mat eigvec;
  bool adjusted=!arma::eig_sym(eigval,eigvec,0.5*(f.h+f.h.t()));
  if(adjusted) { eigval.ones(d); eigvec.eye(d,d); }
  for(int j=0;j<d;++j) if(eigval(j)<curvatureFloor) { eigval(j)=curvatureFloor; adjusted=true; }
  const arma::mat sigma=eigvec*arma::diagmat(1.0/eigval)*eigvec.t();
  double quartic=0.0;
  for(int i=0;i<d;++i) for(int j=0;j<d;++j)
    for(int k=0;k<d;++k) for(int l=0;l<d;++l)
      quartic += f.q.at(f.i4(i,j,k,l))*sigma(i,j)*sigma(k,l);
  arma::mat root;
  arma::chol(root,sigma,"lower");
  std::vector<double> t1(static_cast<std::size_t>(d)*d*d,0.0),
    t2(t1.size(),0.0), tw(t1.size(),0.0);
  auto idx=[d](int i,int j,int k){return (static_cast<std::size_t>(i)*d+j)*d+k;};
  for(int a=0;a<d;++a) for(int j=0;j<d;++j) for(int k=0;k<d;++k)
    for(int i=0;i<d;++i) t1.at(idx(a,j,k))+=root(i,a)*f.t.at(f.i3(i,j,k));
  for(int a=0;a<d;++a) for(int b=0;b<d;++b) for(int k=0;k<d;++k)
    for(int j=0;j<d;++j) t2.at(idx(a,b,k))+=root(j,b)*t1.at(idx(a,j,k));
  for(int a=0;a<d;++a) for(int b=0;b<d;++b) for(int c=0;c<d;++c)
    for(int k=0;k<d;++k) tw.at(idx(a,b,c))+=root(k,c)*t2.at(idx(a,b,k));
  double crossTerm=0.0, traceTerm=0.0;
  for(double value:tw) crossTerm+=value*value;
  for(int a=0;a<d;++a) { double value=0.0; for(int b=0;b<d;++b) value+=tw.at(idx(a,b,b)); traceTerm+=value*value; }
  const double correction=-quartic/8.0+crossTerm/12.0+traceTerm/8.0;
  const double first=-f.v+0.5*d*std::log(2.0*M_PI)-
    0.5*arma::accu(arma::log(eigval));
  const bool reliable=std::isfinite(correction) &&
    std::abs(correction)<=0.5 && 1.0+correction>correctionFloor;
  const double factor=reliable ? 1.0+correction : 1.0;
  if(!reliable) adjusted=true;
  return {first,first+std::log(factor),correction,adjusted,f.g,f.h};
}

} // namespace

// [[Rcpp::export]]
Rcpp::List lmsGraphLaplace2Cpp(const Rcpp::List& matrices,
                               const arma::mat& modes,
                               const int numXis,
                               const int numEtas,
                               const Rcpp::List& dataR,
                               const Rcpp::List& colidxR,
                               const Rcpp::IntegerVector& orderedIndex,
                               const bool logistic=true,
                               const double curvatureFloor=1e-6,
                               const double correctionFloor=1e-8,
                               const int ncores=1) {
  const std::vector<Observation> observations=observationsFromR(dataR,colidxR);
  if(modes.n_rows!=observations.size() || modes.n_cols!=static_cast<arma::uword>(numXis+numEtas))
    Rcpp::stop("Laplace modes have incompatible dimensions.");
  const LaplaceModel model(matrices,numXis,numEtas,orderedIndex,logistic);
  arma::vec first(observations.size()),second(observations.size()),correction(observations.size());
  arma::uvec adjusted(observations.size(),arma::fill::zeros);
  arma::mat modeGradient(observations.size(),numXis+numEtas);
  arma::cube modeHessian(numXis+numEtas,numXis+numEtas,observations.size());
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores>1) schedule(dynamic)
  #endif
  for(int i=0;i<static_cast<int>(observations.size());++i) {
    LaplaceResult value=observationLaplace(model,modes.row(i),observations[i],curvatureFloor,correctionFloor);
    first(i)=value.first; second(i)=value.second; correction(i)=value.correction; adjusted(i)=value.adjusted;
    modeGradient.row(i)=value.gradient.t();
    modeHessian.slice(i)=value.hessian;
  }
  return Rcpp::List::create(Rcpp::Named("first")=first,Rcpp::Named("second")=second,
    Rcpp::Named("correction")=correction,Rcpp::Named("adjusted")=adjusted,
    Rcpp::Named("modeGradient")=modeGradient,
    Rcpp::Named("modeHessian")=modeHessian);
}

// Evaluate a collection of parameter-perturbed models at fixed posterior modes.
// The variants are packed into one parallel loop to avoid repeated R/C++ calls.
// [[Rcpp::export]]
Rcpp::List lmsGraphLaplaceFixedBatchCpp(
    const Rcpp::List& matrixVariants,
    const arma::mat& modes,
    const int numXis,
    const int numEtas,
    const Rcpp::List& dataR,
    const Rcpp::List& colidxR,
    const Rcpp::IntegerVector& orderedIndex,
    const bool secondOrder=true,
    const bool logistic=true,
    const double curvatureFloor=1e-6,
    const double correctionFloor=1e-8,
    const int ncores=1) {
  const std::vector<Observation> observations=observationsFromR(dataR,colidxR);
  const int variants=matrixVariants.size();
  const int dimension=numXis+numEtas;
  if(variants<1) Rcpp::stop("At least one Laplace parameter variant is required.");
  if(modes.n_rows!=observations.size() ||
     modes.n_cols!=static_cast<arma::uword>(dimension))
    Rcpp::stop("Laplace modes have incompatible dimensions.");
  arma::mat contribution(observations.size(),variants);
  arma::cube modeGradient(observations.size(),dimension,variants);
  arma::umat adjusted(observations.size(),variants,arma::fill::zeros);
  std::vector<LaplaceModel> models;
  models.reserve(variants);
  for(int variant=0;variant<variants;++variant) {
    models.emplace_back(Rcpp::as<Rcpp::List>(matrixVariants[variant]),
                        numXis,numEtas,orderedIndex,logistic);
  }
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(ncores) if(ncores>1) schedule(static)
  #endif
  for(int variant=0;variant<variants;++variant) {
    for(int observation=0;
        observation<static_cast<int>(observations.size());++observation) {
      const LaplaceResult value=observationLaplace(
        models[variant],modes.row(observation),observations[observation],
        curvatureFloor,correctionFloor);
      contribution(observation,variant)=secondOrder ? value.second : value.first;
      modeGradient.slice(variant).row(observation)=value.gradient.t();
      adjusted(observation,variant)=value.adjusted;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("contribution")=contribution,
    Rcpp::Named("modeGradient")=modeGradient,
    Rcpp::Named("adjusted")=adjusted);
}
