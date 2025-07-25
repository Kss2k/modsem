createThetaLabel <- function(labelMatrices, labelMatricesCov,
                              constrExprs, start = NULL) {
  matrices <- c(labelMatrices, labelMatricesCov)
  labels <- lapply(matrices, FUN = function(x) {
    if (NCOL(x) == 0 || NROW(x) == 0) return(NULL)

    select <- apply(x, MARGIN = 2, FUN = function(z)
                    !canBeNumeric(z, includeNA = TRUE))
    as.vector(x[select])
  }) |> unlist() |> unique()

  if (!is.null(constrExprs)) {
    labels <- labels[!labels %in% constrExprs$fixedParams]
  }

  if (is.null(start)) {
    start <- vapply(labels, FUN.VALUE = numeric(1L),
                    FUN = function(x) stats::runif(1))
  }

  names(start) <- labels
  start
}


calcThetaLabel <- function(theta, constrExprs) {
  if (is.null(constrExprs)) return(theta)
  Theta <- c(as.list(theta), constrExprs$thetaFixedParams)

  for (i in seq_along(constrExprs$fixedParams)) {
    Theta[constrExprs$fixedParams[[i]]] <-
      eval(constrExprs$exprs[[i]], envir = Theta)
  }

  unlist(Theta)
}


fillMatricesLabels <- function(matrices, labelMatrices, thetaLabel,
                               constraintExprs = NULL) {
  if (is.null(thetaLabel)) return(matrices)
  labels <- names(thetaLabel)
  purrr::map2(.x = matrices, .y = labelMatrices, .f = function(M, L) {
    lapply(labels, FUN = function(label) {
      M[L == label] <<- thetaLabel[[label]]
    })
    M
  })
}


parTableLabelsToList <- function(parTable) {
  parTable <- parTable[parTable$label != "", ]

  labelList <- list()
  for (i in seq_len(NROW(parTable))) {
    val   <- parTable[i, "est"]
    label <- parTable[i, "label"]

    labelList[[label]] <- val
  }

  labelList
}


getVarsExpr <- function(expr) {
  tokens <- createTokensLine(getCharsLine(expr))
  tokens[vapply(tokens, FUN.VALUE = logical(1L), FUN = is.LavName)] |>
    lapply(as.character) |> unlist()
}


removeConstraintExpressions <- function(parTable) {
  if (!NROW(parTable)) return(parTable)
  parTable[!parTable$op %in% c("==", "<", ">", ":="), ]
}


getVCOV_LabelledParams <- function(vcov, model, theta, method, epsilon = 1e-8) {
  J <- getJacobianLabelledParams(model = model, theta = theta, epsilon = epsilon,
                                 method = method)
  if (is.null(J)) return(vcov)
  J %*% vcov %*% t(J)
}


getJacobianLabelledParams <- function(model, theta, method, epsilon = 1e-8) {
  constrExprs <- model$constrExprs

  g <- getTransformationsTheta(model, theta = theta, method = method)
  J <- matrix(0, nrow = length(g), ncol = length(theta),
              dimnames = list(names(g), names(theta)))

  for (i in seq_along(theta)) {
    theta_i      <- theta
    theta_i[[i]] <- theta[[i]] + epsilon
    g_i <- getTransformationsTheta(model, theta = theta_i, method = method)

    J[, i] <- (g_i - g) / epsilon
  }

  J
}


getTransformationsTheta <- function(model, theta, method) {
  theta <- calcPhiTheta(theta = theta, model = model, method = method)
  if (!model$totalLenThetaLabel) return(theta)
  if (model$lenThetaLabel) {
    thetaLabel <- theta[seq_len(model$lenThetaLabel)]
    theta <- theta[-seq_len(model$lenThetaLabel)]
  } else thetaLabel <- NULL

  thetaLabel <- suppressWarnings(calcThetaLabel(thetaLabel, model$constrExprs))
  c(thetaLabel, theta)
}


constraintsContainUnmatchedLabels <- function(parTable, labels) {
  vapply(seq_len(NROW(parTable)), FUN.VALUE = logical(1L), FUN = function(i) {
    if (!parTable[i, "op"] %in% c("==", "<", ">", ":=")) return(FALSE)
    x <- getVarsExpr(parTable[i, "rhs"])
    !all(x %in% labels)
  })
}
