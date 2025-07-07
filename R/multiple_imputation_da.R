#' Estimate a \code{modsem} model using multiple imputation
#'
#' @param model.syntax \code{lavaan} syntax
#'
#' @param data A dataframe with observed variables used in the model.
#'
#' @param method Method to use:
#' \describe{
#'   \item{\code{"lms"}}{latent model structural equations (not passed to \code{lavaan}).}
#'   \item{\code{"qml"}}{quasi maximum likelihood estimation of latent model structural equations (not passed to \code{lavaan}).}
#' }
#'
#' @param m Number of imputations to perform. More imputations will yield better estimates
#'   but can also be (a lot) slower. 
#'
#' @param verbose Should progress be printed to the console?
#'
#' @param se How should corrected standard errors be computed? Alternatives are:
#' \describe{
#'   \item{\code{"simple"}}{Uncorrected standard errors are only calculated once, in the first imputation.
#'                          The standard errors are thereafter corrected using the distribution of the
#'                          estimated coefficients from the different imputations.}
#'   \item{\code{"full"}}{Uncorrected standard errors are calculated and aggregated for each imputation.
#'                        This can give more accurate results, but can be (a lot) slower.
#'                        The standard errors are thereafter corrected using the distribution of the
#'                        estimated coefficients from the different imputations.}
#' }  
#' 
#' @param ... Arguments passed to \code{\link{modsem}}.
#'
#' @details
#'  \code{modsem_impute} is currently only available for the DA approaches 
#'  (LMS and QML). It performs multiple imputation using \code{Amelia::amelia}
#'  and returns aggregated coefficients from the multiple imputations, along with
#'  corrected standard errors.
#'
#' @examples
#' 
#' m1 <- '
#'   # Outer Model
#'   X =~ x1 + x2 +x3
#'   Y =~ y1 + y2 + y3
#'   Z =~ z1 + z2 + z3
#' 
#'   # Inner model
#'   Y ~ X + Z + X:Z
#' '
#' 
#' oneInt2 <- oneInt
#' 
#' set.seed(123)
#' k <- 200
#' I <- sample(nrow(oneInt2), k, replace = TRUE)
#' J <- sample(ncol(oneInt2), k, replace = TRUE)
#' for (k_i in seq_along(I)) oneInt2[I[k_i], J[k_i]] <- NA
#'
#' \dontrun{
#' est <- modsem_mimpute(m1, oneInt2, m = 25)
#' summary(est)
#' }
#' @export
modsem_mimpute <- function(model.syntax, 
                           data, 
                           method = "lms",
                           m = 25, 
                           verbose = interactive(), 
                           se = c("simple", "full"),
                           ...) {
  stopif(!method %in% c("lms", "qml"), 
         "Multiple imputation is only available for the LMS and QML approaches!")

  modsem_mimput_modsem_da(
    model.syntax = model.syntax,
    data = data,
    method = method,
    m = m,
    verbose = verbose,
    ...
  )
}


modsem_mimput_modsem_da <- function(model.syntax,
                                    data,
                                    method = "lms",
                                    m = 25,
                                    cov.syntax = NULL,
                                    verbose = TRUE,
                                    se = c("simple", "full"),
                                    calc.se = NULL,
                                    ...) {
  se <- tolower(se)
  se <- match.arg(se)

  parTable <- rbind(modsemify(model.syntax), modsemify(cov.syntax))
  ovs      <- getOVs(parTable)
  
  missing.cols <- !ovs %in% colnames(data)
  stopif(any(missing.cols), 
         "Variables missing from data:\n", 
         paste(ovs[missing.cols], sep = "\n"))

  data <- data[, ovs, drop = FALSE]

  if (verbose) cat("Imputing data...\n")
  imputed  <- Amelia::amelia(data, m = m, p2s = 0)

  COEF.ALL  <- vector("list", length = m)
  COEF.FREE <- vector("list", length = m)
  VCOV.ALL  <- vector("list", length = m)
  VCOV.FREE <- vector("list", length = m)
  fits <- vector("list", length = m)

  THETA <- NULL

  for (i in seq_len(m)) {
    printedLines <- utils::capture.output(split = TRUE, {
      if (verbose) printf("Fitting imputation %d/%d...\n", i, m)

      if (is.null(calc.se))
        calc.se_i <- ifelse(se == "simple" && i > 1, yes = FALSE, no = TRUE)
      else
        calc.se_i <- calc.se

      if (is.null(THETA)) {
        optimize <- TRUE
        start    <- NULL
      } else {
        optimize <- FALSE
        start    <- apply(THETA, MARGIN = 2, FUN = mean, na.rm = TRUE) # na.rm should't be necessary, ever...
      }

      data_i <- as.data.frame(imputed$imputations[[i]])

      fit_i <- tryCatch(
        modsem_da(
          model.syntax = model.syntax,
          data         = data_i,
          method       = method,
          cov.syntax   = cov.syntax,
          start        = start,
          verbose      = verbose,
          optimize     = optimize,
          calc.se      = calc.se_i,
          ...
        ), error = \(e) {print(e); NULL}
      )

      if (is.null(fit_i)) next


      fits[[i]]      <- fit_i
      COEF.ALL[[i]]  <- coef(fit_i, type = "all")
      COEF.FREE[[i]] <- coef(fit_i, type = "free")

      if (i > 1 && se == "simple") {
        VCOV.ALL[[i]]  <- VCOV.ALL[[1]]
        VCOV.FREE[[i]] <- VCOV.FREE[[1]]

      } else if (calc.se_i) {
        VCOV.ALL[[i]]  <- vcov(fit_i, type = "all")
        VCOV.FREE[[i]] <- vcov(fit_i, type = "free")

      } else {
        k.all  <- length(COEF.ALL[[i]])
        k.free <- length(COEF.FREE[[i]])
        d.all  <- names(COEF.ALL[[i]])
        d.free <- names(COEF.FREE[[i]])

        VCOV.ALL[[i]]  <- matrix(0, nrow = k.all, ncol = k.all, 
                                 dimnames = list(d.all, d.all))
        VCOV.FREE[[i]] <- matrix(0, nrow = k.free, ncol = k.free, 
                                 dimnames = list(d.free, d.free))
      }

      THETA <- rbind(
        THETA, 
        matrix(fit_i$theta, nrow = 1, dimnames = list(NULL, names(fit_i$theta)))
      )
    })

    nprinted <- length(printedLines)
    if (i < m) eraseConsoleLines(nprinted)
  }

  failed <- vapply(fits, FUN.VALUE = logical(1L), FUN = is.null)
 
  if (any(failed)) {
    warning2(sprintf("Model estimation failed in %d out of %d impuations!", 
                     sum(failed), m), immediate. = FALSE)

    fits      <- fits[!failed]
    COEF.ALL  <- COEF.ALL[!failed]
    COEF.FREE <- COEF.FREE[!failed]
    VCOV.ALL  <- VCOV.ALL[!failed]
    VCOV.FREE <- VCOV.FREE[!failed]
    m         <- sum(!failed)
  }

  pool.all  <- rubinPool(COEF.ALL,  VCOV.ALL)
  pool.free <- rubinPool(COEF.FREE, VCOV.FREE)

  coef.all  <- pool.all$theta.bar
  vcov.all  <- pool.all$Tvcov
  
  coef.free <- pool.free$theta.bar
  vcov.free <- pool.free$Tvcov

  parTable1   <- parameter_estimates(fits[[1]])
  orig.labels <- parTable1$label
  parTable1   <- getMissingLabels(parTable1)
  parTableT   <- data.frame(label = names(coef.all), est.t = coef.all)

  parTable    <- leftJoin(x = parTable1, y = parTableT, by = "label")
  match       <- !is.na(parTable$est.t)

  parTable$est[match] <- parTable$est.t[match]
  parTable$est.t      <- NULL
  parTable$label[!parTable$label %in% orig.labels] <- ""

  matrices    <- aggregateMatrices(fits, type = "main")
  covMatrices <- aggregateMatrices(fits, type = "cov")
  expected.matrices <- aggregateMatrices(fits, type = "expected")

  getScalarFit <- function(fit, field) 
    vapply(fits, FUN.VALUE = numeric(1L), \(fit) fit[[field]])

  fit.out <- fits[[1]]
  fit.out$coefs.all      <- coef.all
  fit.out$coefs.free     <- coef.free
  fit.out$vcov.all       <- vcov.all
  fit.out$vcov.free      <- vcov.free 
  fit.out$parTable       <- parTable
  fit.out$information    <- sprintf("Rubin-corrected (m=%d)", m)
  fit.out$FIM            <- solve(vcov.free)
  fit.out$theta          <- apply(THETA, MARGIN = 2, FUN = mean, na.rm = TRUE)
  fit.out$iterations     <- sum(getScalarFit(fits, field = "iterations"))
  fit.out$logLik         <- mean(getScalarFit(fits, field = "logLik"))

  fit.out$model$matrices          <- matrices
  fit.out$model$covModel$matrices <- covMatrices
  fit.out$expected.matrices       <- expected.matrices

  fit.out$imputations <- list(fitted = fits, data = imputed)

  fit.out
}


matrixMeans <- function(matrices, na.rm = TRUE) {
  if (!length(matrices)) 
    return(matrix(nrow=0, ncol=0))

  # Take the mean of a list of n*p matrices
  # With the option of ignoring missing values

  flattened <- as.data.frame(lapply(matrices, FUN = c))
  means <- rowMeans(flattened, na.rm = na.rm)

  dims <- dimnames(matrices[[1]])
  n    <- NROW(matrices[[1]])
  p    <- NCOL(matrices[[1]])

  matrix(means, nrow = n, ncol = p, dimnames = dims)
}


aggregateMatrices <- function(fits, type) {
  switch(type,
    main = {
      matrices  <- fits[[1]]$model$matrices
      names.num <- namesParMatrices
    },
    cov = {
      matrices  <- fits[[1]]$model$covModel$matrices
      names.num <- namesParMatricesCov
    },
    expected = {
      matrices  <- fits[[1]]$expected.matrices
      names.num <- names(matrices)
    }
  )

  if (is.null(matrices) || !length(matrices))
    return(matrices)

  m <- length(fits)

  for (i in seq_len(m)[-1]) {
    fit_i      <- fits[[i]]
    matrices_i <- switch(type,
      main     = fit_i$model$matrices,
      cov      = fit_i$model$covModel$matrices,
      expected = fit_i$expected.matrices
    )
    matrices_i <- matrices_i[names(matrices_i) %in% names.num]

    for (mat in names(matrices_i)) {
      if (is.null(matrices_i[[mat]]) || !length(matrices_i[[mat]]))
        next

      matrices[[mat]] <- matrices[[mat]] + matrices_i[[mat]] 
    }
  }

  for (mat in names.num)
    matrices[[mat]] <- matrices[[mat]] / m

  matrices
}


rubinPool <- function(coef.list, vcov.list) {
  stopifnot(length(coef.list) == length(vcov.list))

  m     <- length(coef.list)
  THETA <- do.call(cbind, coef.list)
  p     <- nrow(THETA)

  # Pooled point estimate
  theta.bar <- rowMeans(THETA, na.rm = TRUE) # na.rm shouldn't be necessary
  # but just in case...

  # Within-imputation covariance (W)
  W <- matrixMeans(vcov.list, na.rm = TRUE)

  # Between-imputation covariance (B)
  THETAc <- sweep(THETA, 1, theta.bar, "-")
  B      <- (THETAc %*% t(THETAc)) / (m - 1)

  # Total Rubin-corrected covariance (T)
  Tvcov  <- W + (1 + 1/m) * B

  # Small-sample df
  df <- (m - 1) * (1 + diag(W) / ((1 + 1/m) * diag(B)))^(-2)

  list(theta.bar = theta.bar,
       W         = W,
       B         = B,
       Tvcov     = Tvcov,
       df        = df)
}
