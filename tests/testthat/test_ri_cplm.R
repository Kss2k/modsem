devtools::load_all()

model.inp <- '
  # Create between components (random intercepts)
  RIemo =~ 1 * emo_3 + 1 * emo_5 + 1 * emo_7;
  RIcon =~ 1 * con_3 + 1 * con_5 + 1 * con_7;


  # Create within-person centered variables
  wemo_3 =~ 1 * emo_3;
  wemo_5 =~ 1 * emo_5;
  wemo_7 =~ 1 * emo_7;

  wcon_3 =~ 1 * con_3;
  wcon_5 =~ 1 * con_5;
  wcon_7 =~ 1 * con_7;
                    
  # moderator also needs to be decomposed into within and between parts 	              
  RIper =~ 1 * per_3 + 1 * per_5 + 1 * per_7;
            
  wper_3 =~ 1 * per_3;
  wper_5 =~ 1 * per_5;
  wper_7 =~ 1 * per_7;
                    
  # Constrain the measurement error variances to zero
 	# or close to zero to allow for reasobale imputation times
  con_3 ~~ 0.2 * con_3
  con_5 ~~ 0.2 * con_5
  con_7 ~~ 0.2 * con_7
  emo_3 ~~ 0.2 * emo_3
  emo_5 ~~ 0.2 * emo_5
  emo_7 ~~ 0.2 * emo_7
  per_3 ~~ 0.2 * per_3
  per_5 ~~ 0.2 * per_5
  per_7 ~~ 0.2 * per_7
	            
  # Estimate the covariance between the random intercepts
  RIemo ~~ RIper
  RIemo ~~ RIcon
  RIper ~~ RIcon

  # Estimate the lagged effects between
  # the within-person centered variables
  wemo_7 ~ wemo_5 + wper_5 + wcon_5
  wper_7 ~ wemo_5 + wper_5 + wcon_5
  wcon_7 ~ wemo_5 + wper_5 + wcon_5

  wemo_5 ~ wemo_3 + wper_3 + wcon_3
  wper_5 ~ wemo_3 + wper_3 + wcon_3
  wcon_5 ~ wemo_3 + wper_3 + wcon_3
                    
  # Specify interaction terms between within-person centred 
	# conduct problems and the within-person centered peer problems
  # predict within-person centred emotional problems with interaction
  wemo_7 ~ wcon_5:wper_5
  wemo_5 ~ wcon_3:wper_3

  # Estimate the covariance between the within-person
  # components at the first wave
  wemo_3 ~~ wper_3
  wemo_3 ~~ wcon_3
  wper_3 ~~ wcon_3

  # Estimate the covariances between the residuals of
  # the within-person components (the innovations)
  wemo_5 ~~ wper_5
  wemo_5 ~~ wcon_5
  wper_5 ~~ wcon_5

  wemo_7 ~~ wper_7
  wemo_7 ~~ wcon_7
  wper_7 ~~ wcon_7
'

model.out <- '
  # Create between components (random intercepts)
  RIemo =~ 1 * emo_3 + 1 * emo_5 + 1 * emo_7;
  RIcon =~ 1 * con_3 + 1 * con_5 + 1 * con_7;

  # Create within-person centered variables
  wemo_3 =~ 1 * emo_3;
  wemo_5 =~ 1 * emo_5;
  wemo_7 =~ 1 * emo_7;

  wcon_3 =~ 1 * con_3;
  wcon_5 =~ 1 * con_5;
  wcon_7 =~ 1 * con_7;
                    
  # moderator also needs to be decomposed into within and between parts 	              
  RIper =~ 1 * per_3 + 1 * per_5 + 1 * per_7;
            
  wper_3 =~ 1 * per_3;
  wper_5 =~ 1 * per_5;
  wper_7 =~ 1 * per_7;
	            
  wemo_7 ~ 0.361 * wemo_5 + 0.097 * wper_5 + 0.147 * wcon_5 + 0.105 * wcon_5:wper_5
  wper_7 ~ 0.030 * wemo_5 + 0.348 * wper_5 + 0.123 * wcon_5
  wcon_7 ~ 0.074 * wemo_5 + 0.056 * wper_5 + 0.086 * wcon_5

  wemo_5 ~ 0.142 * wemo_3 + 0.071 * wper_3 + 0.086 * wcon_3 + (-0.010) * wcon_3:wper_3
  wper_5 ~ 0.018 * wemo_3 + 0.100 * wper_3 + 0.050 * wcon_3
  wcon_5 ~ (-0.008) * wemo_3 + 0.006 * wper_3 + 0.105 * wcon_3

  RIemo ~~  0.481 * RIper + 0.438 * RIcon
  RIcon ~~  0.469 * RIper
  wemo_3 ~~ 0.328 * wper_3 + 0.478 * wcon_3
  wper_3 ~~ 0.427 * wcon_3
  wemo_5 ~~ 0.397 * wper_5 + 0.301 * wcon_5
  wper_5 ~~ 0.202 * wcon_5
  wemo_7 ~~ 0.521 * wper_7
  wemo_7 ~~ 0.475 * wcon_7
  wper_7 ~~ 0.311 * wcon_7

  con_3 ~ 2.827 * 1
  con_5 ~ 1.517 * 1
  con_7 ~ 1.402 * 1
  emo_3 ~ 1.385 * 1
  emo_5 ~ 1.407 * 1
  emo_7 ~ 1.528 * 1
  per_3 ~ 1.563 * 1
  per_5 ~ 1.176 * 1
  per_7 ~ 1.242 * 1

  RIemo  ~~ 0.806 * RIemo 
  RIper  ~~ 0.758 * RIper 
  RIcon  ~~ 1.290 * RIcon 
  wemo_3 ~~ 1.234 * wemo_3
  wper_3 ~~ 1.614 * wper_3
  wcon_3 ~~ 2.743 * wcon_3

  con_3  ~~ 0.200 * con_3  
  con_5  ~~ 0.200 * con_5  
  con_7  ~~ 0.200 * con_7  
  emo_3  ~~ 0.200 * emo_3  
  emo_5  ~~ 0.200 * emo_5  
  emo_7  ~~ 0.200 * emo_7  
  per_3  ~~ 0.200 * per_3  
  per_5  ~~ 0.200 * per_5  
  per_7  ~~ 0.200 * per_7  
  wemo_5 ~~ 1.489 * wemo_5 
  wemo_7 ~~ 1.808 * wemo_7 
  wper_5 ~~ 1.143 * wper_5 
  wper_7 ~~ 1.287 * wper_7 
  wcon_5 ~~ 0.764 * wcon_5 
  wcon_7 ~~ 0.881 * wcon_7 
'

# Simulate data
parTableSim <- modsemify(model.out)
parTableSim$est <- as.numeric(parTableSim$mod)

set.seed(238947)
data <- modsem:::simulateDataParTable(
  parTable = parTableSim,
  N = 10000
)$OV[[1L]]

fit.lms <- modsem(
  model.syntax = model.inp,
  data         = data,
  method       = "lms",
  nodes        = 15,
  adaptive.quad = "aghq",
  auto.split.syntax = TRUE,
  optimize = FALSE,
  verbose  = TRUE
)
summary(fit.lms, H0 = FALSE)
