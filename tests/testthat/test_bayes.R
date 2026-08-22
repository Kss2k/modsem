devtools::load_all()

model <- '
  X =~ x1 + b * x2 + x3
  Z =~ z1 + b * z2 + x3 # let x3 cross-load
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
  b := normal(0.0, 1.2)
'

parTable <- modsemify(model, parentheses.as.string = TRUE)
parTableEx <- expandModsemParTable(parTable)

cat(buildStanSyntaxFromParTable(parTableEx))


m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

parTable <- modsemify(m1)
parTableEx <- expandModsemParTable(parTable)

library(rstan)
rstan_options(auto_write = TRUE, threads_per_chain = 4)      # cache compiled models
options(mc.cores = parallel::detectCores()) 

stanCode <- buildStanSyntaxFromParTable(parTableEx)
stanModel <- rstan::stan_model(model_code = stanCode)

stanData <- list(
  N = nrow(oneInt),
  IND__x1 = oneInt$x1,
  IND__x2 = oneInt$x2,
  IND__x3 = oneInt$x3,
  IND__z1 = oneInt$z1,
  IND__z2 = oneInt$z2,
  IND__z3 = oneInt$z3,
  IND__y1 = oneInt$y1,
  IND__y2 = oneInt$y2,
  IND__y3 = oneInt$y3
)

fit <- rstan::sampling(
  object  = stanModel,
  data    = stanData,
  chains  = 2,
  iter    = 4000,
  warmup  = 2000,
  pars    = c("MAT__ZETA"),
  include = FALSE
)data {
int<lower=0> N;
vector[N] IND__x1;
vector[N] IND__x2;
vector[N] IND__x3;
vector[N] IND__z1;
vector[N] IND__z2;
vector[N] IND__z3;
vector[N] IND__y1;
vector[N] IND__y2;
vector[N] IND__y3;
}

parameters {
vector[N] MAT__ZETA[3];
real X__INTR__1;
real Z__INTR__1;
real Y__INTR__1;
real X__COV__X;
real Z__COV__Z;
real Y__COV__Y;
real Z__COV__X;
real Y__REG__X;
real Y__REG__Z;
real Y__REG__X__MOD__Z;
real IND__x1__INTR__1;
real IND__x1__COV__IND__x1;
real IND__x2__INTR__1;
real IND__x2__COV__IND__x2;
real X__MSR__x2;
real IND__x3__INTR__1;
real IND__x3__COV__IND__x3;
real X__MSR__x3;
real IND__z1__INTR__1;
real IND__z1__COV__IND__z1;
real IND__z2__INTR__1;
real IND__z2__COV__IND__z2;
real Z__MSR__z2;
real IND__z3__INTR__1;
real IND__z3__COV__IND__z3;
real Z__MSR__z3;
real IND__y1__INTR__1;
real IND__y1__COV__IND__y1;
real IND__y2__INTR__1;
real IND__y2__COV__IND__y2;
real Y__MSR__y2;
real IND__y3__INTR__1;
real IND__y3__COV__IND__y3;
real Y__MSR__y3;
}

transformed parameters {
matrix[3,3] MAT__PSI;
vector[3] VEC_ZETA_MU;
VEC_ZETA_MU[1] = X__INTR__1;
VEC_ZETA_MU[2] = Z__INTR__1;
VEC_ZETA_MU[3] = Y__INTR__1;
MAT__PSI[1,1] = X__COV__X;
MAT__PSI[1,1] = X__COV__X;
MAT__PSI[2,2] = Z__COV__Z;
MAT__PSI[2,2] = Z__COV__Z;
MAT__PSI[3,3] = Y__COV__Y;
MAT__PSI[3,3] = Y__COV__Y;
MAT__PSI[2,1] = Z__COV__X;
MAT__PSI[1,2] = Z__COV__X;
vector[N] LV__X = MAT__ZETA[1];
vector[N] LV__Z = MAT__ZETA[2];
vector[N] RES__Y = MAT__ZETA[3];
vector[N] LV__Y = RES__Y + Y__REG__X*LV__X + Y__REG__Z*LV__Z + Y__REG__X__MOD__Z*LV__X.*LV__Z;
real X__MSR__x1 = 1;
vector[N] RES__IND__x1 = IND__x1 - (IND__x1__INTR__1 + X__MSR__x1*LV__X);
vector[N] RES__IND__x2 = IND__x2 - (IND__x2__INTR__1 + X__MSR__x2*LV__X);
vector[N] RES__IND__x3 = IND__x3 - (IND__x3__INTR__1 + X__MSR__x3*LV__X);
real Z__MSR__z1 = 1;
vector[N] RES__IND__z1 = IND__z1 - (IND__z1__INTR__1 + Z__MSR__z1*LV__Z);
vector[N] RES__IND__z2 = IND__z2 - (IND__z2__INTR__1 + Z__MSR__z2*LV__Z);
vector[N] RES__IND__z3 = IND__z3 - (IND__z3__INTR__1 + Z__MSR__z3*LV__Z);
real Y__MSR__y1 = 1;
vector[N] RES__IND__y1 = IND__y1 - (IND__y1__INTR__1 + Y__MSR__y1*LV__Y);
vector[N] RES__IND__y2 = IND__y2 - (IND__y2__INTR__1 + Y__MSR__y2*LV__Y);
vector[N] RES__IND__y3 = IND__y3 - (IND__y3__INTR__1 + Y__MSR__y3*LV__Y);
}

model {
MAT__ZETA ~ multi_normal(VEC_ZETA_MU, MAT__PSI);
RES__IND__x1 ~ normal(0.0, sqrt(IND__x1__COV__IND__x1));
RES__IND__x2 ~ normal(0.0, sqrt(IND__x2__COV__IND__x2));
RES__IND__x3 ~ normal(0.0, sqrt(IND__x3__COV__IND__x3));
RES__IND__z1 ~ normal(0.0, sqrt(IND__z1__COV__IND__z1));
RES__IND__z2 ~ normal(0.0, sqrt(IND__z2__COV__IND__z2));
RES__IND__z3 ~ normal(0.0, sqrt(IND__z3__COV__IND__z3));
RES__IND__y1 ~ normal(0.0, sqrt(IND__y1__COV__IND__y1));
RES__IND__y2 ~ normal(0.0, sqrt(IND__y2__COV__IND__y2));
RES__IND__y3 ~ normal(0.0, sqrt(IND__y3__COV__IND__y3));
}

