# Global variables
namesParMatrices <- c("lambdaX", "lambdaY", "gammaXi", "gammaEta",
                      "thetaDelta", "thetaEpsilon", "phi", "A",
                      "psi", "tauX", "tauY", "alpha", "beta0", "omegaEtaXi",
                      "omegaXiXi")
namesParMatricesCov <- c("gammaXi", "gammaEta", "A", "psi", "phi")


createTheta <- function(model, start = NULL, parTable.in = NULL) {
  THETA_LAB  <- NULL
  THETA_COV  <- NULL
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
    thetaLabel_g <- createThetaLabel(submodel$labelMatrices,
                                     submodel$covModel$labelMatrices,
                                     model$params$constrExprs)
    THETA_LAB <- unionByNames(THETA_LAB, thetaLabel_g)
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
    labelsMain_g <- LABELS_MAIN_GROUPS[[g]]
    labelsCov_g  <- LABELS_COV_GROUPS[[g]]

    selectTL <- seq_along(THETA_LAB) # available to all sub models
    selectTC <- which(names(THETA_COV) %in% labelsCov_g)
    selectTM <- which(names(THETA_MAIN) %in% labelsMain_g)
    
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


createThetaCovModel <- function(covModel, start = NULL) {
  M <- covModel$matrices

  phi      <- as.vector(M$phi)
  psi      <- as.vector(M$psi)
  alpha    <- as.vector(M$alpha)
  gammaXi  <- as.vector(M$gammaXi)
  gammaEta <- as.vector(M$gammaEta)
  thetaCov <- c("phi" = phi,
                "psi" = psi,
                "gammaXi" = gammaXi,
                "gammaEta" = gammaEta)

  lavLabelsCov <- createLavLabelsCov(M, subset = is.na(thetaCov))
  thetaCov <- thetaCov[is.na(thetaCov)]
  thetaCov <- fillThetaIfStartNULL(start = start, theta = thetaCov,
                                   lavlab = lavLabelsCov)

  list(theta = thetaCov, lavLabels = lavLabelsCov)
}


fillThetaIfStartNULL <- function(start,
                                 theta,
                                 lavlab = NULL,
                                 var = 1,
                                 cov = 0,
                                 meas = 0.7,
                                 mean = 0,
                                 reg = 0) {
  if (!is.null(start)) {
    return(theta)

  } else if (!is.null(lavlab)) {
    tryCatch({
      OP <- "~~|=~|~1|~"
      op <- stringr::str_extract(lavlab, pattern = OP)
      lr <- stringr::str_split_fixed(lavlab, pattern = OP, n = 2)

      lhs <- lr[, 1]
      rhs <- lr[, 2]
      op[is.na(op)] <- "~"

      theta.filled                          <- theta
      theta.filled[op == "~"]               <- reg
      theta.filled[op == "=~"]              <- meas
      theta.filled[op == "~1"]              <- mean
      theta.filled[op == "~~" & lhs == rhs] <- var
      theta.filled[op == "~~" & lhs != rhs] <- cov
      theta.filled[is.na(theta.filled)]     <- reg

      theta.filled
    }, error = \(e)
        fillThetaIfStartNULL(start = start, theta = theta, lavlab = NULL)
    )
  } else {
    vapply(theta, FUN = function(x) stats::runif(1),
           FUN.VALUE = numeric(1L))
  }
}


fillModel <- function(model, theta, fillPhi = FALSE, method = "lms") {
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

    # cov model
    thetaCov <- NULL
    if (length(params.utils$SELECT_THETA_COV[[g]]))
      thetaCov <- theta[params.utils$SELECT_THETA_COV[[g]]]

    thetaMain <- theta[params.utils$SELECT_THETA_MAIN[[g]]]

    submodel$covModel <- fillCovModel(submodel$covModel, thetaCov, thetaLabel)
    submodel$matrices <- fillMainModel(submodel, thetaMain, thetaLabel,
                                       fillPhi = fillPhi, method = method)

    model$models[[g]] <- submodel
  }

  model
}


fillMainModel <- function(model, theta, thetaLabel, fillPhi = FALSE,
                          method = "lms") {
  xis      <- model$info$xis
  numXis   <- model$info$numXis
  numEtas  <- model$info$numEtas
  M        <- model$matrices
  covModel <- model$covModel

  lMatrices <- model$labelMatrices[namesParMatrices]
  pMatrices <- M[namesParMatrices]
  M[namesParMatrices] <- fillMatricesLabels(pMatrices, lMatrices, thetaLabel)

  if (!is.null(model$covModel$matrices)) {
    M$phi <- M$A <- expectedCovModel(covModel, method = method, sortedXis = xis)
  } else if (method == "lms") {
    M$A <- fillNA_Matrix(M$A, theta = theta, pattern = "^A([0-9]*)")
  } else if (method == "qml") {
    M$phi <- fillSymmetric(M$phi, fetch(theta, "^phi"))
  }

  M$lambdaX      <- fillNA_Matrix(M$lambdaX, theta = theta, pattern = "^lambdaX")
  M$lambdaY      <- fillNA_Matrix(M$lambdaY, theta = theta, pattern = "^lambdaY")
  M$thetaDelta   <- fillSymmetric(M$thetaDelta, fetch(theta, "^thetaDelta"))
  M$thetaEpsilon <- fillSymmetric(M$thetaEpsilon, fetch(theta, "thetaEpsilon"))
  M$psi          <- fillSymmetric(M$psi, fetch(theta, "^psi"))
  M$tauX         <- fillNA_Matrix(M$tauX, theta = theta, pattern = "^tauX")
  M$tauY         <- fillNA_Matrix(M$tauY, theta = theta, pattern = "^tauY")
  M$alpha        <- fillNA_Matrix(M$alpha, theta = theta, pattern = "^alpha")
  M$beta0        <- fillNA_Matrix(M$beta0, theta = theta, pattern = "^beta0")
  M$gammaEta     <- fillNA_Matrix(M$gammaEta, theta = theta, pattern = "^gammaEta")
  M$gammaXi      <- fillNA_Matrix(M$gammaXi, theta = theta, pattern = "^gammaXi")
  M$omegaXiXi    <- fillNA_Matrix(M$omegaXiXi, theta = theta, pattern = "^omegaXiXi")
  M$omegaEtaXi   <- fillNA_Matrix(M$omegaEtaXi, theta = theta, pattern = "^omegaEtaXi")

  if (fillPhi) M$phi <- M$A %*% t(M$A)
  M
}


fillCovModel <- function(covModel, theta, thetaLabel) {
  if (is.null(names(theta))) names(theta) <- names(covModel$theta)
  if (is.null(covModel$matrices)) return(covModel)
  M <- covModel$matrices

  lMatrices <- covModel$labelMatrices[namesParMatricesCov]
  pMatrices <- M[namesParMatricesCov]
  M[namesParMatricesCov] <- fillMatricesLabels(pMatrices, lMatrices, thetaLabel)

  M$psi      <- fillSymmetric(M$psi, fetch(theta, "^psi"))
  M$gammaEta <- fillNA_Matrix(M$gammaEta, theta = theta, pattern = "^gammaEta")
  M$gammaXi  <- fillNA_Matrix(M$gammaXi, theta = theta, pattern = "^gammaXi")
  M$phi <- fillSymmetric(M$phi, fetch(theta, "^phi"))

  covModel$matrices <- M
  covModel
}


fillNA_Matrix <- function(X, theta, pattern) {
  idx <- is.na(X) & !is.nan(X)
  values <- fetch(theta, pattern)
  if (length(values) && sum(idx) != length(values)) {
    stop2("Mismatch when filling matrix for pattern `", pattern, "`: expected ",
          sum(idx), " values but got ", length(values))
  }
  if (sum(idx) > 0 && length(values) == 0) {
    stop2("No values found in theta vector for pattern `", pattern, "`.")
  }
  X[idx] <- values
  X
}


fillSymmetric <- function(mat, values) {
  idx <- is.na(mat) & !is.nan(mat)
  if (sum(idx) > 0 && length(values) == 0) {
    stop2("No values provided to fill symmetric matrix.")
  }
  if (length(values) && length(values) != sum(idx)) {
    stop2("Mismatch when filling symmetric matrix: expected ", sum(idx),
          " values but got ", length(values))
  }
  mat[idx] <- values
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  mat
}


getParamBounds <- function(model) {
  lower <- rep(-Inf, model$params$freeParams)
  upper <- rep(Inf, model$params$freeParams)
  names(lower) <- names(upper) <- names(model$theta)

  parTable <- model$parTable
  BOUNDS <- parTable[canBeNumeric(parTable$rhs) & parTable$op %in%
                     BOUNDUARY_OPS, , drop = FALSE]
  bound.param <- BOUNDS$lhs
  bound.type  <- BOUNDS$op
  bound.value <- as.numeric(BOUNDS$rhs)

  upper.type    <- bound.type == "<"
  lower.type <- bound.type == ">"

  bound.upper <- structure(bound.value[upper.type],
                           names = bound.param[upper.type])
  bound.lower <- structure(bound.value[lower.type],
                           names = bound.param[lower.type])

  upper[names(bound.upper)] <- bound.upper
  lower[names(bound.lower)] <- bound.lower

  list(lower = lower, upper = upper)
}


checkStartingParams <- function(start, model) {
  if (length(start) != length(model$theta)) {
    stop2("The length of the starting parameters does not match the number of parameters in the model")
  }
  if (is.null(names(start))) {
    names(start) <- names(model$theta)
  }
  if (!all(names(start) %in% names(model$theta))) {
    stop2("The names of the starting parameters do not match the names of the parameters in the model")
  }

  NULL
}


calcPhiTheta <- function(theta, model, method) {
  if (method != "lms") return(theta)

  modFilled <- fillModel(theta = theta, model = model, method = method,
                         fillPhi = TRUE)

  for (g in seq_len(model$info$n.groups)) {
    select <- c(model$params$SELECT_THETA_LAB[[g]],
                model$params$SELECT_THETA_COV[[g]],
                model$params$SELECT_THETA_MAIN[[g]])

    theta_g      <- theta[select]
    submodFilled <- modFilled$models[[g]]
    submodel     <- model$models[[g]]

    if (!is.null(submodel$covModel$matrices)) {
      matEst <- submodFilled$covModel$matrices
      matLab <- submodel$covModel$labelMatrices
      matNA  <- submodel$covModel$matrices
    } else {
      matEst <- submodFilled$matrices
      matNA  <- submodel$matrices
      matLab <- submodel$labelMatrices
    }

    vals   <- as.vector(matEst$phi[is.na(matNA$A) & !is.nan(matNA$A)])
    labels <- as.vector(matLab$A)

    if (any(labels != "")) {
      allVals <- as.vector(matEst$phi)
      labVals <- allVals[labels != ""]
      labels  <- labels[labels != ""]
      missing <- setdiff(labels, names(theta))
      stopif(length(missing),
             "Missing labelled parameters in theta vector: ",
             paste(missing, collapse = ", "))
      
      theta[select][labels] <- labVals
    }

    theta[select][grepl("^A([0-9]*)", names(theta[select]))] <- vals
  }

  theta
}


LMS_BLOCKS = list(
  lambdaX      = 0,
  lambdaY      = 1,
  tauX         = 2,
  tauY         = 3,
  thetaDelta   = 4,
  thetaEpsilon = 5,
  A            = 6,
  psi          = 7,
  alpha        = 8,
  beta0        = 9,
  gammaXi      = 10,
  gammaEta     = 11,
  omegaXiXi    = 12,
  omegaEtaXi   = 13,
  phi          = NA
)


SYMMETRIC_BLOCKS_LMS = c(
  thetaDelta = 4,
  thetaEpsilon = 5,
  psi = 7,
  phi = NA
)


getParamNamesMatrix <- function(mat, matname) {
  if (is.character(mat)) {
    c(mat)

  } else {
    M <- list()
    M[[matname]] <- mat
    names(unlist(M))
  }
}


getParamLocationsMatrices <- function(matrices, isFree = is.na, g = 1L, ignore.g.label = FALSE) {
  matrices <- matrices[intersect(names(matrices), names(LMS_BLOCKS))]
  locations <- data.frame(param = NULL, block = NULL, row = NULL, col = NULL)
  for (blockname in names(matrices)) {
    X <- matrices[[blockname]]
    n <- nrow(X)
    m <- ncol(X)

    if (!any(isFree(X))) next

    params <- getParamNamesMatrix(mat = X, matname = blockname)
    block  <- LMS_BLOCKS[[blockname]]
    rowidx <- matrix(seq_len(n) - 1, nrow = n, ncol = m, byrow = FALSE)
    colidx <- matrix(seq_len(m) - 1, nrow = n, ncol = m, byrow = TRUE)

    params <- params[isFree(X)]
    rowidx <- rowidx[isFree(X)]
    colidx <- colidx[isFree(X)]

    if (g > 1L && !ignore.g.label)
      params <- sprintf("%s.g%d", params, g)

    locationsBlock <- data.frame(
      param = params,
      group = g,
      block = block,
      row   = rowidx,
      col   = colidx,
      symmetric = as.integer(block %in% SYMMETRIC_BLOCKS_LMS)
    )

    locations <- rbind(locations, locationsBlock)
  }

  locations
}


getGradientStruct <- function(model, theta) {
  tryCatch(
    getGradientStructSimple(model = model, theta = theta),
    error = function(e) {
      warning2("Failed to compute gradient structure: ", e$message)

      list(
        locations   = NULL,
        Jacobian    = NULL,
        nlinDerivs  = NULL,
        evalTheta   = NULL,
        hasCovModel = TRUE, # may not be true, but we should behave as if it is
        isNonLinear = TRUE  # may not be true, but we should behave as if it is
      )
    }
  )
}


getGradientStructSimple <- function(model, theta) {
  hasCovModel <- !is.null(model$models[[1L]]$covModel$matrices)

  if (hasCovModel) {
    out <- list(
      locations   = NULL,
      Jacobian    = NULL,
      nlinDerivs  = NULL,
      evalTheta   = NULL,
      hasCovModel = TRUE,
      isNonLinear = TRUE  # may not be true, but we should behave as if it is
    )

    return(out)
  }

  parTable <- model$parTable
  parTable <- parTable[!parTable$op %in% BOUNDUARY_OPS, , drop = FALSE] # not relevant

  customParams <- parTable[parTable$op %in% c(":=", "=="), , drop = FALSE]
  for (i in seq_len(NROW(customParams))) {
    row <- customParams[i, , drop = FALSE]
    eq  <- sprintf("(%s)", row$rhs)
    pattern <- sprintf("(?<![A-z_\\.])%s(?![A-z_\\.])", row$lhs)

    mask <- parTable$op %in% CONSTRAINT_OPS
    parTable$rhs[mask] <- stringr::str_replace_all(
      parTable$rhs[mask], pattern = pattern, replacement = eq
    )
  }

  isConstraint <- parTable$op %in% CONSTRAINT_OPS & !canBeNumeric(parTable$rhs)
  constraints  <- parTable[isConstraint, ]
  restParTable <- parTable[!isConstraint, ]
  # TODO: find better way of removing unecessary constraints
  constraints  <- constraints[constraints$lhs %in% restParTable$mod, ]

  derivatives <- list()
  derivatives2 <- list()
  for (i in seq_len(NROW(constraints))) {
    constrVar <- constraints[i, "lhs"]
    constrEq  <- constraints[i, "rhs"]

    derivatives[[constrVar]] <- derivateConstraint(constrEq)
    derivatives2[[constrVar]] <- secondDerivateConstraint(constrEq)
  }

  isLinear <- vapply(
    X = derivatives, FUN.VALUE = logical(1L),
    FUN = \(X) all(vapply(X, FUN.VALUE = logical(1L), FUN = is.atomic))
  )

  linDerivs   <- derivatives[isLinear]
  nlinDerivs  <- derivatives[!isLinear]
  nlinDerivs2 <- derivatives2[!isLinear]
  evalTheta   <- \(theta) c(theta, suppressWarnings(calcThetaLabel(theta, model$constrExprs))) # This could be made a bit better

  locations <- NULL
  for (g in seq_len(model$info$n.groups)) {
    submodel <- model$models[[g]]

    locations <- rbind(
       locations,
       getParamLocationsMatrices(submodel$matrices, isFree = is.na, g = g),
       getParamLocationsMatrices(submodel$labelMatrices, isFree = \(x) x != "",
                                 g = g, ignore.g.label = TRUE) # labels don't change with g here
    )
  }

  k <- nrow(locations)
  m <- length(theta)

  param.full <- locations$param
  param.part <- names(theta)

  ordering <- structure(seq_along(theta), names = param.part)
  ordering <- ordering[param.full]

  locations  <- locations[order(ordering), ]
  param.full <- locations$param

  Jacobian <- matrix(0, nrow = m, ncol = k,
                     dimnames = list(param.part, param.full))
  Jacobian2 <- Jacobian

  for (par in param.full) {
    match.full <- param.full == par
    match.part <- param.part == par

    Jacobian[match.part, match.full] <- 1
  }

  for (dep in names(linDerivs)) {
    deriv <- linDerivs[[dep]]

    for (indep in names(deriv)) {
      match.full <- param.full == dep
      match.part <- param.part == indep
      Jacobian[match.part, match.full] <- deriv[[indep]]
    }
  }

  # In multigroup models we can have duplicated labels, that doesn't work
  # with out logic, so we must add some unique identifiers to the param.full
  # We also pass param.full as is, to map back and forth.
  addid <- \(nm) paste0(nm, "#", seq_along(nm))
  colnames(Jacobian)  <- addid(colnames(Jacobian))
  colnames(Jacobian2) <- addid(colnames(Jacobian2))

  locations$param <- addid(locations$param)

  list(
    param.full  = param.full,
    locations   = locations,
    Jacobian    = Jacobian,
    Jacobian2   = Jacobian2,
    nlinDerivs  = nlinDerivs,
    nlinDerivs2 = nlinDerivs2,
    evalTheta   = evalTheta,
    hasCovModel = hasCovModel,
    isNonLinear = length(nlinDerivs) > 1
  )
}


derivateConstraint <- function(constr) {
  f <- stats::formula(paste0("~", constr))
  eq <- Deriv::Deriv(f, combine = "list", cache.exp = FALSE)

  vars <- all.vars(f)
  k    <- length(vars)

  if (k == 0) return(NULL)

  if (is.null(names(eq)) && k == 1L) {
    eq <- stats::setNames(list(eq), nm = vars)
  } else if (is.null(names(eq))) {
    names(eq) <- all.vars(f)
  } else {
    eq <- as.list(eq)[-1]
  }

  eq
}


secondDerivateConstraint <- function(constr) {
  f <- stats::formula(paste0("~", constr))
  eq <- Deriv::Deriv(f, nderiv = 2, combine = "list", cache.exp = FALSE)

  vars <- all.vars(f)
  k    <- length(vars)

  if (k == 0) return(NULL)

  out <- list()
  if (is.null(names(eq)) && k == 1L) {
    eq <- stats::setNames(list(eq), nm = vars)
    return(eq)

  } else if (is.null(names(eq))) {
    names(eq) <- all.vars(f)
    return(eq)

  } else {
    for (indep in all.vars(f)) {
      out[[indep]] <- eq[[indep]][[indep]]
    }
    return(out)
  }
}
