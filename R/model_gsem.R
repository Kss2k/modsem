createGsemModelGroup <- function(parTable, ordered = NULL) {
  # General Information
  higherOrderLVs <- getHigherOrderLVs(parTable)
  indsHigherOrderLVs <- getIndsLVs(parTable, lVs = higherOrderLVs, isOV = FALSE)
  ovs <- getOVs(parTable)

  # endogenous variables (etas)
  etas    <- getSortedEtas(parTable, isLV = TRUE, checkAny = TRUE)
  numEtas <- length(etas)

  indsEtas    <- getIndsLVs(parTable, lVs = etas, isOV = TRUE, ovs = ovs)
  numIndsEtas <- vapply(indsEtas, FUN.VALUE = vector("integer", 1L),
                        FUN = length)
  allIndsEtas    <- unique(unlist(indsEtas))
  numAllIndsEtas <- length(allIndsEtas)

  # exogenouts variables (xis) and interaction terms
  intTerms      <- getIntTermRows(parTable)
  varsInts      <- getVarsInts(intTerms)
  allVarsInInts <- unique(unlist(varsInts))
  xis           <- getXis(parTable, checkAny = TRUE)
  numXis        <- length(xis)

  indsXis    <- getIndsLVs(parTable, lVs = xis, isOV = TRUE, ovs = ovs)
  numIndsXis <- vapply(indsXis, FUN.VALUE = vector("integer", 1L),
                       FUN = length)
  allIndsXis    <- unique(unlist(indsXis))
  numAllIndsXis <- length(allIndsXis)

  xwith <- unique(stringr::str_replace_all(intTerms$rhs, ":", "__XWITH__"))
  lvs <- c(xwith, xis, etas)
  indsLVs <- c(indsXis, indsEtas)
  numLVs <- length(lvs)

  # clean data
  data.cleaned <- prepDataModsemDA(data, allIndsXis, allIndsEtas,
                                   missing = missing, cluster = cluster)

  OMEGA <- matrix(0, nrow = xwith, ncol = numEtas + numXis,
                  dimnames = list(xwith, c(xis, etas)))
  for (i in seq_along(xwith)) {
    vars <- stringr::str_split(intTerms$rhs[[i]], pattern = ":")[[1L]]
    OMEGA[i, colnames(OMEGA) %in% vars] <- 1L
  }

  listLambda <- constructLambda(lvs, indsLVs,
                                parTable = parTable,
                                auto.fix.first = auto.fix.first)

  lambda      <- listLambda$numeric
  labelLambda <- listLambda$label

  listTau <- constructTau(lvs, indsLVs,
                          parTable = parTable, mean.observed = mean.observed)

  tau      <- listTau$numeric
  labelTau <- listTau$label

  listTheta <- constructTheta(lvs, indsLVs, parTable = parTable,
                              auto.fix.single = auto.fix.single)

  theta      <- listTheta$numeric
  thetaLabel <- listTheta$label


  # structural model
  Ieta         <- diag(numLVs) # used for (B^-1 = (Ieta - gammaEta)^-1)
  listGammaXi  <- constructGamma(lvs, lvs, parTable = parTable,
                                 auto.fix.first = auto.fix.first)
  gamma      <- listGamma$numeric
  labelGamma <- listGamma$label

  # covariance matrices
  listPsi  <- constructPsi(lvs, parTable = parTable, orthogonal.y = orthogonal.y)
  psi      <- listPsi$numeric
  labelPsi <- listPsi$label

  listPhiXi <- constructPhi(xis, method = method, cov.syntax = cov.syntax,
                            parTableCovModel = parTableCovModel,
                            parTable = parTable, orthogonal.x = orthogonal.x)
  phiXi      <- listPhiXi$numeric
  labelPhiXi <- listPhiXi$label

  psi[xis, xis]      <- phiXi
  labelPsi[xis, xis] <- labelPhiXi
  psi[xwith, xwith]  <- diag(length(xwith))
  labelPsi[xwith, xwith] <- ""

  # mean etas
  listAlpha <- constructAlpha(lvs, parTable = parTable,
                              mean.observed = mean.observed)
  alpha      <- listAlpha$numeric
  labelAlpha <- listAlpha$label

  listThresholds  <- constructThresholds(ordered = ordered, parTable = parTable, data = data)
  thresholds      <- listThresholds$numeric
  labelThresholds <- listThresholds$label

  is.ordered <- rownames(lambda) %in% ordered

  # list of matrices
  matrices <- list(
    lambdaX      = lambda,
    gammaEta     = gamma,
    thetaEpsilon = theta,
    Ieta         = Ieta,
    psi          = psi,
    tauX         = tau,
    alpha        = alpha,
  )

  labelMatrices <- list(
    lambdaX      = lambda,
    gammaEta     = gamma,
    thetaEpsilon = theta,
    Ieta         = Ieta,
    psi          = psi,
    tauX         = tau,
    alpha        = alpha,
  )

  quad <- quadrature(m, k = numEtas + numXis, quad.range = quad.range,
                     adaptive = adaptive.quad, adaptive.frequency = adaptive.frequency)

  model <- list(
    info = list(
      N             = data.cleaned$n,
      ncol          = data.cleaned$k,
      xis           = xis,
      etas          = etas,
      numXis        = numXis,
      numEtas       = numEtas,
      indsXis       = indsXis,
      indsEtas      = indsEtas,
      allIndsXis    = allIndsXis,
      allIndsEtas   = allIndsEtas,
      varsInts      = varsInts,
      mean.observed = mean.observed,
      is.ordered    = is.ordered,

      has.interaction    = NROW(intTerms) > 0L,
      higherOrderLVs     = higherOrderLVs,
      indsHigherOrderLVs = indsHigherOrderLVs,

      lavOptimizerSyntaxAdditions = lavOptimizerSyntaxAdditions
    ),

    data          = data.cleaned,
    data.raw      = data,
    quad          = quad,
    matrices      = matrices,
    labelMatrices = labelMatrices,
    syntax        = syntax,
    cov.syntax    = cov.syntax,
    parTable      = parTable,
    covModel      = covModel,
    lavaan.fit    = NULL
  )

  if (checkModel)
    preCheckModel(model = model, covModel = covModel, method = method,
                  missing = missing)

  model
}


specifyModelGsem <- function(group.info, ...) {


}
# Global variables
namesParMatricesGsem <- c("lambda", "gamma",
                          "theta", "phi", "tau", "alpha")


createThetaGsem <- function(model, start = NULL, parTable.in = NULL) {
  THETA_LAB  <- NULL
  THETA_MAIN <- NULL

  n.groups   <- model$info$n.groups
  if (is.null(n.groups)) n.groups <- model$info$n.groups
  if (is.null(n.groups)) n.groups <- length(model$models)

  LABELS_MAIN_GROUPS <- vector("list", length = n.groups)
  LABELS_COV_GROUPS  <- vector("list", length = n.groups)

  unionByNames <- function(x, y) {
    new <- setdiff(names(y), names(x))
    if (length(new)) c(x, y[new]) else x
  }

  # The custom labels must all be added, before computing THETA_LAB_ALL
  for (g in seq_len(n.groups)) {
    submodel     <- model$models[[g]]
    thetaLabel.g <- createThetaLabel(submodel$labelMatrices,
                                     submodel$covModel$labelMatrices,
                                     model$params$constrExprs)
    THETA_LAB <- unionByNames(THETA_LAB, thetaLabel.g)
  }
  THETA_LAB_ALL <- calcThetaLabel(THETA_LAB, model$params$constrExprs)
  LAV_LAB       <- names(THETA_LAB_ALL)

  for (g in seq_len(n.groups)) {
    submodel <- model$models[[g]]
    etas <- submodel$info$etas

    listThetaCov <- createThetaCovModel(submodel$covModel)
    thetaCov     <- listThetaCov$theta
    lavLabelsCov <- listThetaCov$lavLabels

    M            <- submodel$matrices
    lambdaX      <- as.vector(M$lambdaX)
    lambdaY      <- as.vector(M$lambdaY)
    thetaDelta   <- as.vector(M$thetaDelta)
    thetaEpsilon <- as.vector(M$thetaEpsilon)
    phi          <- as.vector(M$phi)
    A            <- as.vector(M$A)
    psi          <- as.vector(M$psi)
    tauX         <- as.vector(M$tauX)
    tauY         <- as.vector(M$tauY)
    alpha        <- as.vector(M$alpha)
    beta0        <- as.vector(M$beta0)
    gammaXi      <- as.vector(M$gammaXi)
    gammaEta     <- as.vector(M$gammaEta)
    omegaXiXi    <- as.vector(M$omegaXiXi)
    omegaEtaXi   <- as.vector(M$omegaEtaXi)

    allModelValues <- c("lambdaX"      = lambdaX,
                        "lambdaY"      = lambdaY,
                        "tauX"         = tauX,
                        "tauY"         = tauY,
                        "thetaDelta"   = thetaDelta,
                        "thetaEpsilon" = thetaEpsilon,
                        "phi"          = phi,
                        "A"            = A,
                        "psi"          = psi,
                        "alpha"        = alpha,
                        "beta0"        = beta0,
                        "gammaXi"      = gammaXi,
                        "gammaEta"     = gammaEta,
                        "omegaXiXi"    = omegaXiXi,
                        "omegaEtaXi"   = omegaEtaXi)

    lavLabelsMain <- createLavLabels(M, subset = is.na(allModelValues),
                                     etas = etas, parTable.in = parTable.in)

    thetaMain <- allModelValues[is.na(allModelValues)]
    thetaMain <- fillThetaIfStartNULL(start = start, theta = thetaMain,
                                      lavlab = lavLabelsMain)

    allLabels <- names(c(thetaCov, thetaMain))
    lavLabels <- combineLavLabels(lavLabelsMain = lavLabelsMain,
                                  lavLabelsCov = lavLabelsCov,
                                  currentLabels = allLabels,
                                  g = g)

    if (g > 1L) {
      .newnames <- \(nm) sprintf("%s.g%d", nm, g)

      # thetaLabel is labelled across the submodels, so the names don't change!
      names(thetaCov)    <- .newnames(names(thetaCov))
      names(thetaMain)   <- .newnames(names(thetaMain))
    }

    LABELS_MAIN_GROUPS[[g]] <- if (length(thetaMain)) names(thetaMain) else character(0)
    LABELS_COV_GROUPS[[g]]  <- if (length(thetaCov)) names(thetaCov) else character(0)

    THETA_COV  <- unionByNames(THETA_COV, thetaCov)
    THETA_MAIN <- unionByNames(THETA_MAIN, thetaMain)
    LAV_LAB    <- union(LAV_LAB, lavLabels)
  }

  THETA <- c(THETA_LAB, THETA_COV, THETA_MAIN)

  SELECT_THETA_LAB  <- vector("list", length = n.groups)
  SELECT_THETA_COV  <- vector("list", length = n.groups)
  SELECT_THETA_MAIN <- vector("list", length = n.groups)

  for (g in seq_len(n.groups)) {
    labelsMain.g <- LABELS_MAIN_GROUPS[[g]]
    labelsCov.g  <- LABELS_COV_GROUPS[[g]]

    selectTL <- seq_along(THETA_LAB) # available to all sub models
    selectTC <- which(names(THETA_COV) %in% labelsCov.g)
    selectTM <- which(names(THETA_MAIN) %in% labelsMain.g)

    # The selections must be offset by their respective location THETA as a whole
    # only selectTL doesn't need to be shifted, as it's at the start
    selectTC <- selectTC + length(THETA_LAB)
    selectTM <- selectTM + length(THETA_LAB) + length(THETA_COV)

    SELECT_THETA_LAB[[g]]  <- selectTL
    SELECT_THETA_COV[[g]]  <- selectTC
    SELECT_THETA_MAIN[[g]] <- selectTM
  }

  list(theta             = THETA,
       freeParams        = length(THETA),
       SELECT_THETA_LAB  = SELECT_THETA_LAB,
       SELECT_THETA_COV  = SELECT_THETA_COV,
       SELECT_THETA_MAIN = SELECT_THETA_MAIN,
       lavLabels         = LAV_LAB)
}


specifyModelGsem <- function(..., group.info, createTheta = TRUE) {
  args <- list(...)

  n.groups <- group.info$n.groups
  stopif(n.groups < 1L, "Invalid grouping structure supplied.")

  parTable    <- group.info$parTable
  parTableCov <- group.info$parTableCov
  group.col   <- parTable$group

  stopif(is.null(group.col) || max(group.col) != n.groups,
         "Number of group-specific parameter tables does not match number of groups.")

  submodels <- vector("list", length = n.groups)

  for (g in seq_len(n.groups)) {
    args.g <- args

    if (!is.null(group.info$data))
      args.g$data <- group.info$data[group.info$indices[[g]], , drop = FALSE]
    else
      args.g <- addNamedNullField(args.g, field = "data")

    args.g$parTable <- parTable[parTable$group == g, , drop = FALSE]

    if (!is.null(parTableCov))
      args.g$parTableCovModel <- parTableCov[parTableCov$group == g, , drop = FALSE]
    else
      args.g <- addNamedNullField(args.g, field = "parTableCovModel")

    submodel.g <- do.call(specifModelDA_Group, args.g)
    submodel.g$info$group <- group.info$levels[[g]]
    submodel.g$info$n.groups <- 1L

    submodels[[g]] <- submodel.g
  }

  model <- list(
    models   = submodels,
    syntax   = args$syntax,
    data.raw = group.info$data.raw,
    parTable = parTable,
    info     = list(
      n.groups      = n.groups,
      group.levels  = group.info$levels,
      group.info    = group.info,

      # Constants across groups
      xis           = submodels[[1L]]$info$xis,
      etas          = submodels[[1L]]$info$etas,
      numXis        = submodels[[1L]]$info$numXis,
      numEtas       = submodels[[1L]]$info$numEtas,
      indsXis       = submodels[[1L]]$info$indsXis,
      indsEtas      = submodels[[1L]]$info$indsEtas,
      allIndsXis    = submodels[[1L]]$info$allIndsXis,
      allIndsEtas   = submodels[[1L]]$info$allIndsEtas,
      varsInts      = submodels[[1L]]$info$varsInts,
      latentEtas    = submodels[[1L]]$info$latentEtas,
      scalingInds   = submodels[[1L]]$info$scalingInds,
      kOmegaEta     = submodels[[1L]]$info$kOmegaEta,
      nonLinearXis  = submodels[[1L]]$info$nonLinearXis,
      mean.observed = submodels[[1L]]$info$mean.observed,

      has.interaction    = submodels[[1L]]$info$has.interaction,
      higherOrderLVs     = submodels[[1L]]$info$higherOrderLVs,
      indsHigherOrderLVs = submodels[[1L]]$info$indsHigherOrderLVs,

      lavOptimizerSyntaxAdditions = submodels[[1L]]$lavOptimizerSyntaxAdditions
    ),

    params = list()
  )

  # Currenlty we assume covModel has an uniform structure
  model$params$constrExprs <- getConstrExprs(parTable, model$models[[1L]]$covModel$parTable)

  if (createTheta) {
    params <- createTheta(model, parTable.in = parTable)
    model$params[names(params)] <- params

    # TODO: Remove occurences of `model$theta`, and replace them with `model$params$theta`
    model$theta <- params$theta # an ugly design decision, that was made at the very start

    model$params$bounds <- getParamBounds(model)
    model$params$gradientStruct <- getGradientStruct(model, theta = params$theta)
  }

  model
}


fillModelGsem <- function(model, theta, fillPhi = FALSE, method = "lms") {
  params.utils <- model$params

  if (is.null(names(theta)))
    names(theta) <- names(params.utils$theta)

  # labeled parameters
  thetaLabel <- NULL

  if (length(params.utils$SELECT_THETA_LAB)) {
    thetaLabel <- theta[params.utils$SELECT_THETA_LAB[[1L]]] # same for all groups
    thetaLabel <- suppressWarnings(calcThetaLabel(thetaLabel, params.utils$constrExprs))
  }

  for (g in seq_len(model$info$n.groups)) {
    submodel <- model$models[[g]]

    thetaMain <- theta[params.utils$SELECT_THETA_MAIN[[g]]]
    submodel$matrices <- fillGroupModelGsem(submodel, thetaMain, thetaLabel,
                                            fillPhi = fillPhi, method = method)

    model$models[[g]] <- submodel
  }

  model
}


fillGroupModelGsem <- function(model, theta, thetaLabel) ) {
  M        <- model$matrices

  lMatrices <- model$labelMatrices[namesParMatricesGsem]
  pMatrices <- M[namesParMatrices]
  M[namesParMatrices] <- fillMatricesLabels(pMatrices, lMatrices, thetaLabel)

  M$lambdaX      <- fillNA_Matrix(M$lambdaX, theta = theta, pattern = "^lambdaX")
  M$thetaDelta   <- fillSymmetric(M$thetaDelta, fetch(theta, "^thetaDelta"))
  M$psi          <- fillSymmetric(M$psi, fetch(theta, "^psi"))
  M$tauX         <- fillNA_Matrix(M$tauX, theta = theta, pattern = "^tauX")
  M$alpha        <- fillNA_Matrix(M$alpha, theta = theta, pattern = "^alpha")
  M$gammaXi      <- fillNA_Matrix(M$gammaXi, theta = theta, pattern = "^gammaXi")

  M
}
