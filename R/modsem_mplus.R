#' Estimation latent interactions through \code{Mplus}
#'
#' @param model.syntax lavaan/modsem syntax. Multilevel models can be specified
#'   using \code{level: 1} (or \code{level: within}) and \code{level: 2}
#'   (or \code{level: between}) blocks, see the details and examples.
#'
#' @param data dataset
#'
#' @param estimator estimator argument passed to \code{Mplus}. Bayesian estimation is
#'   available using \code{estimator = "bayes"}, see the details.
#'
#' @param cluster cluster argument passed to \code{Mplus}. Required for multilevel
#'   models. For single level models it is used to compute standard errors
#'   corrected for clustering (i.e., \code{TYPE = COMPLEX}).
#'
#' @param categorical categorical argument passed to \code{Mplus}.
#'
#' @param type type argument passed to \code{Mplus}. If \code{NULL} (default) it is
#'   chosen based on the model: \code{"twolevel"} for multilevel models, \code{"complex"}
#'   if \code{cluster} is specified, and \code{"general"} otherwise. \code{"random"} is
#'   added if the model contains interaction effects or random slopes.
#'
#' @param algorithm algorithm argument passed to \code{Mplus}. If \code{NULL} (default)
#'   it is set to \code{"integration"} for models with interaction effects (unless
#'   \code{estimator = "bayes"}), otherwise it is left to \code{Mplus}.
#'
#' @param processors processors argument passed to \code{Mplus}.
#'
#' @param integration integration argument passed to \code{Mplus}. Only used if
#'   \code{algorithm = "integration"}.
#'
#' @param analysis.options Additional options passed to the \code{ANALYSIS} command of
#'   \code{Mplus}, as a named list (or character vector). E.g.,
#'   \code{list(BITERATIONS = "10000 (2000)", CHAINS = 4, BSEED = 123)} for
#'   \code{estimator = "bayes"}.
#'
#' @param rcs Should latent variable indicators be replaced with reliability-corrected
#'   single item indicators instead? See \code{\link{relcorr_single_item}}.
#'   Not available for multilevel models.
#'
#' @param rcs.choose Which latent variables should get their indicators replaced with
#'   reliability-corrected single items? It is passed to \code{\link{relcorr_single_item}}
#'   as the \code{choose} argument.
#'
#' @param rcs.scale.corrected Should reliability-corrected items be scale-corrected? If \code{TRUE}
#'   reliability-corrected single items are corrected for differences in factor loadings between
#'   the items. Default is \code{TRUE}.
#'
#' @param output.std Should \code{STANDARDIZED} be added to \code{OUTPUT}? It is
#'   ignored if \code{Mplus} doesn't support standardized estimates for the model
#'   (e.g., \code{TYPE = TWOLEVEL RANDOM} with a frequentist estimator).
#'
#' @param within Variables which should be declared as \code{WITHIN} variables in
#'   multilevel models. If \code{NULL} (default) all observed variables which only
#'   appear in the within (\code{level: 1}) model are declared as \code{WITHIN} variables.
#'
#' @param between Variables which should be declared as \code{BETWEEN} variables in
#'   multilevel models. If \code{NULL} (default) all observed variables which only
#'   appear in the between (\code{level: 2}) model are declared as \code{BETWEEN} variables.
#'
#' @param cleanup Should the \code{Mplus} files (\code{.inp}, \code{.dat} and \code{.out})
#'   created during estimation be deleted when \code{modsem_mplus()} exits? Default is
#'   \code{FALSE}, meaning the files are kept in the working directory.
#'
#' @param ... arguments passed to other functions
#'
#' @details
#' Interaction terms (e.g., \code{X:Z}) involving at least one latent variable are
#' specified using \code{XWITH}, whilst interaction terms between observed variables are
#' computed as product terms in the \code{DEFINE} command of \code{Mplus}.
#'
#' If \code{estimator = "bayes"} the parameter table contains the posterior means
#' (\code{est}), posterior standard deviations (\code{std.error}), the one-tailed
#' p-values reported by \code{Mplus} (\code{p.value}), and the 95\% credible
#' intervals (\code{ci.lower} and \code{ci.upper}). \code{z.value} is \code{NA}.
#' The number of iterations, chains etc. can be set using \code{analysis.options}.
#'
#' In multilevel models random slopes are specified using the \code{rv()} modifier
#' in the within level model (e.g., \code{Y ~ rv("s")*X}), which is translated to
#' \code{s | Y ON X;} in \code{Mplus}. The random slope (\code{s}) is a latent variable
#' at the between level, where it can be used like any other latent variable. Note that
#' the fixed effect of a random slope is the mean (or intercept) of the random slope at
#' the between level.
#'
#' @return modsem_mplus object
#' @export
#'
#' @examples
#' # Theory Of Planned Behavior
#' m1 <- '
#' # Outer Model
#'   X =~ x1 + x2
#'   Z =~ z1 + z2
#'   Y =~ y1 + y2
#'
#' # Inner model
#'   Y ~ X + Z + X:Z
#' '
#'
#' \dontrun{
#' # Check if Mplus is installed
#' run <- tryCatch({MplusAutomation::detectMplus(); TRUE},
#'                 error = \(e) FALSE)
#'
#' if (run) {
#'   est_mplus <- modsem_mplus(m1, data = oneInt)
#'   summary(est_mplus)
#' }
#'
#' # Multilevel model with a random slope, and a (within level) interaction
#' m2 <- '
#'   level: 1
#'     fw =~ y1 + y2 + y3
#'     fw ~ rv("s")*x1 + x2 + x3 + x2:x3
#'   level: 2
#'     fb =~ y1 + y2 + y3
#'     fb + s ~ w1 + w2
#' '
#'
#' if (run) {
#'   est_mplus_ml <- modsem_mplus(m2, data = lavaan::Demo.twolevel, cluster = "cluster")
#'   summary(est_mplus_ml)
#' }
#' }
#'
modsem_mplus <- function(model.syntax,
                         data,
                         estimator = "ml",
                         cluster = NULL,
                         type = NULL,
                         algorithm = NULL,
                         processors = 2,
                         integration = 15,
                         analysis.options = NULL,
                         rcs = FALSE,
                         rcs.choose = NULL,
                         rcs.scale.corrected = TRUE,
                         output.std = TRUE,
                         categorical = NULL,
                         within = NULL,
                         between = NULL,
                         cleanup = FALSE,
                         ...) {
  if (!requireNamespace("MplusAutomation", quietly = TRUE))
    mod_msg_stop("modsem_mplus() requires the 'MplusAutomation' package; please install it.")

  if (rcs) { # use reliability-correct single items?
    corrected <- relcorr_single_item(
      syntax          = model.syntax,
      data            = data,
      choose          = rcs.choose,
      scale.corrected = rcs.scale.corrected,
      warn.lav        = FALSE
    )

    model.syntax <- corrected$syntax
    data         <- corrected$data
  }

  spec <- specifyMplusModel(
    parTable    = modsemify(model.syntax),
    data        = data,
    cluster     = cluster,
    categorical = categorical,
    within      = within,
    between     = between
  )

  analysis <- specifyMplusAnalysis(
    spec        = spec,
    estimator   = estimator,
    type        = type,
    algorithm   = algorithm,
    processors  = processors,
    integration = integration,
    options     = analysis.options
  )

  input <- MplusAutomation::mplusObject(
    TITLE           = "Running Model via Mplus",
    VARIABLE        = spec$VARIABLE,
    DEFINE          = spec$DEFINE,
    ANALYSIS        = analysis$ANALYSIS,
    MODEL           = spec$MODEL,
    MODELCONSTRAINT = spec$MODELCONSTRAINT,
    OUTPUT          = specifyMplusOutput(spec, analysis = analysis, output.std = output.std),
    usevariables    = colnames(spec$data),
    rdata           = spec$data
  )

  model <- runMplusModel(input, cleanup = cleanup)

  parTable <- mplusParTable(model$results$parameters$unstandardized, spec = spec)
  parTable <- addMplusParamIds(parTable, TECH1 = model$results$tech1, spec = spec)
  free     <- getFreeParamsMplus(parTable)
  coefs    <- structure(parTable$est[free$row], names = free$label)
  vcov     <- getVcovMplus(mode, freee = free)

  modelSpec <- list(
    parTable = modsemParTable(removeParamIdsMplus(parTable)),
    model    = model,
    data     = data,
    coefs    = modsemVector(coefs),
    vcov     = modsemMatrix(vcov, symmetric = TRUE),
    info     = list(
      spec      = spec,
      analysis  = analysis,
      cluster   = cluster,
      intTerms  = spec$intTerms
    )
  )

  structure(modelSpec,
            class = "modsem_mplus",
            method = "Mplus")
}


MPLUS_CONSTRAINT_OPS <- c(":=", "==", "<", ">")
MPLUS_MAX_NAME_LENGTH <- 8L
MPLUS_MAX_LINE_WIDTH <- 85L # The limit is 90 characters, but we keep a margin
MPLUS_RESERVED_WORDS <- c(
  "BY", "ON", "WITH", "PON", "PWITH", "XWITH", "IS", "ARE", "ALL", "TO", "AT",
  "AND", "OR", "NOT", "EQ", "NE", "GT", "LT", "GE", "LE"
)


# ──────────────────────────────────────────────────────────────────────────────
# Model Specification
# ──────────────────────────────────────────────────────────────────────────────
specifyMplusModel <- function(parTable,
                              data,
                              cluster = NULL,
                              categorical = NULL,
                              within = NULL,
                              between = NULL) {
  data <- as.data.frame(data)

  isConstraint <- parTable$op %in% MPLUS_CONSTRAINT_OPS
  constraints  <- parTable[isConstraint, , drop = FALSE]
  parTable     <- addLevelsParTable(parTable[!isConstraint, , drop = FALSE])

  multilevel <- !all(parTable$level == "")
  levels     <- unique(parTable$level)

  mod_stopif(multilevel && is.null(cluster),
    "The `cluster` argument must be specified for multilevel models!"
  )

  mod_stopif(!is.null(cluster) && !all(cluster %in% colnames(data)),
    sprintf("Cluster variable `%s` not found in data!", cluster)
  )

  slopes <- getRandomSlopesMplus(parTable, multilevel = multilevel)

  # Variables by level
  # NOTE: For single level models the level is "", so lists by level are
  # indexed by position, rather than by name
  vars <- lapply(levels, FUN = \(level) {
    getVariablesLevelMplus(
      parTable = parTable[parTable$level == level, , drop = FALSE],
      slopes   = slopes$name
    )
  })

  observed <- unique(unlist(lapply(vars, \(v) c(v$observed, v$elements))))
  missing  <- setdiff(observed, colnames(data))
  mod_stopif(length(missing),
    sprintf("Missing variables in data: %s", paste0(missing, collapse = ", "))
  )

  # Mplus names are case insensitive, can contain at most 8 characters,
  # and can only contain letters, numbers and underscores.
  labels <- getLabelsMplus(parTable, constraints = constraints)
  names  <- getMplusNames(c(
    observed, cluster,
    unlist(lapply(vars, \(v) v$latent)),
    slopes$name, labels)
  )

  intTerms <- vector("list", length(levels))
  for (i in seq_along(levels)) {
    intTerms[[i]] <- specifyInteractionsMplus(vars[[i]], names = names)
    names <- c(names, intTerms[[i]]$names)
  }

  # Observed variables used in the model at each level (including product terms)
  used <- lapply(seq_along(levels), FUN = function(i) {
    unique(c(vars[[i]]$observed, intTerms[[i]]$observed, intTerms[[i]]$products))
  })
  names(used) <- levels

  # Product terms are computed by Mplus (DEFINE), and are not part of the data
  products <- unique(unlist(lapply(intTerms, \(x) x$products)))
  defined  <- unique(unlist(lapply(intTerms, \(x) x$DEFINE)))

  # Data
  isFactor    <- vapply(data, FUN.VALUE = logical(1L), FUN = is.factor)
  categorical <- union(categorical, colnames(data)[isFactor])
  categorical <- intersect(categorical, unlist(lapply(vars, \(v) v$endogenous)))

  data <- data[c(observed, cluster)]
  data <- lapplyDf(data, FUN = \(x) if (is.factor(x)) as.integer(x) else x)

  usevars <- setdiff(unlist(used), products)
  data     <- data[unique(c(usevars, cluster))]
  colnames(data) <- unname(names[colnames(data)])

  # MODEL
  MODEL <- vapply(seq_along(levels), FUN.VALUE = character(1L), FUN = function(i) {
    header <- c(within = "%WITHIN%", between = "%BETWEEN%")[levels[[i]]]

    statements <- c(
      if (!is.na(header)) header,
      intTerms[[i]]$XWITH,
      parTableToMplusModel(
        parTable[parTable$level == levels[[i]], , drop = FALSE],
        names = names
      )
    )

    paste0(statements, collapse = "\n")
  })

  # VARIABLE
  levelVars <- getLevelVariablesMplus(
    used = used,
    within = within,
    between = between
  )

  # Variables created using DEFINE must be listed last in USEVARIABLES
  VARIABLE <- c(
    if (length(products)) mplusVariableList("USEVARIABLES", names[c(usevars, products)]),
    if (!is.null(cluster)) sprintf("CLUSTER = %s;", names[[cluster]]),
    mplusVariableList("CATEGORICAL", names[categorical]),
    mplusVariableList("WITHIN", names[levelVars$within]),
    mplusVariableList("BETWEEN", names[levelVars$between])
  )

  list(
    data            = data,
    names           = names,
    multilevel      = multilevel,
    levels          = levels,
    vars            = vars,
    intTerms        = unique(unlist(lapply(vars, \(v) v$intTerms))),
    slopes          = slopes,
    parTable        = parTable,
    constraints     = constraints,
    has.cluster     = !is.null(cluster),
    has.xwith       = any(vapply(intTerms, \(x) length(x$XWITH) > 0L, logical(1L))),
    VARIABLE        = paste0(VARIABLE, collapse = "\n"),
    DEFINE          = if (length(defined)) paste0(defined, collapse = "\n"),
    MODEL           = paste0(MODEL, collapse = "\n"),
    MODELCONSTRAINT = parTableToMplusConstraints(constraints, names = names)
  )
}


addLevelsParTable <- function(parTable) {
  isLevel <- parTable$lhs == "level" & parTable$op == ":"

  if (!any(isLevel)) {
    parTable$level <- rep("", NROW(parTable))
    return(parTable)
  }

  mod_stopif(!isLevel[[1L]],
    "Multilevel models must start with a `level:` block!"
  )

  level  <- vapply(parTable$rhs[isLevel], FUN.VALUE = character(1L), FUN = function(lvl) {
    switch(tolower(lvl),
      "1"       = "within",
      "2"       = "between",
      "within"  = "within",
      "between" = "between",
       mod_msg_stop(sprintf(
          "Unrecognized level `%s`, expected `1` (within) or `2` (between)!", lvl
       ))
    )
  })

  mod_stopif(anyDuplicated(level),
    "Each level can only be specified once in multilevel models!"
  )

  parTable$level <- unname(level)[cumsum(isLevel)]
  parTable[!isLevel, , drop = FALSE]
}


getRandomSlopesMplus <- function(parTable, multilevel) {
  pattern  <- "^rv\\(\\s*[\"']?([^\"')]+)[\"']?\\s*\\)$"
  isRandom <- grepl(pattern, parTable$mod)

  slopes <- parTable[isRandom, c("lhs", "op", "rhs", "level"), drop = FALSE]
  slopes$name <- stringr::str_match(parTable$mod[isRandom], pattern)[, 2L]

  mod_stopif(NROW(slopes) && !multilevel,
    "Random slopes (`rv()`) are only available in multilevel models!"
  )
  mod_stopif(any(slopes$level != "within"),
    "Random slopes (`rv()`) must be specified in the within level (`level: 1`)!"
  )
  mod_stopif(any(!slopes$op %in% c("~", "=~")),
    "Random slopes (`rv()`) can only be specified for `~` and `=~` parameters!"
  )
  mod_stopif(anyDuplicated(slopes$name),
    "Each random slope (`rv()`) must have a unique name!"
  )

  slopes
}


getVariablesLevelMplus <- function(parTable, slopes = NULL) {
  intTerms <- unique(c(parTable$lhs[grepl(":", parTable$lhs)],
                       parTable$rhs[grepl(":", parTable$rhs)]))
  latent   <- unique(parTable$lhs[parTable$op == "=~"])

  # The random slopes are latent variables at the between level (and not
  # variables at the within level)
  latent <- c(latent, intersect(slopes, c(parTable$lhs, parTable$rhs)))

  hasRhsVar <- !parTable$op %in% c("~1", "|")
  variables <- unique(c(parTable$lhs, parTable$rhs[hasRhsVar]))
  elements  <- unique(unlist(stringr::str_split(intTerms, ":")))

  isIndicator <- parTable$op == "=~"
  isOutcome   <- parTable$op %in% c("~", "|")
  endogenous  <- c(parTable$rhs[isIndicator], parTable$lhs[isOutcome])

  list(
    intTerms   = intTerms,
    latent     = unique(latent),
    observed   = setdiff(variables, c(latent, intTerms)),
    elements   = setdiff(elements, latent), # observed variables in interaction terms
    endogenous = setdiff(endogenous, c(latent, intTerms))
  )
}


getLabelsMplus <- function(parTable, constraints) {
  isLabel <- parTable$mod != "" & !canBeNumeric(parTable$mod) &
    !grepl("^rv\\(", parTable$mod) & parTable$mod != "NA"

  unique(c(parTable$mod[isLabel], constraints$lhs[constraints$op == ":="]))
}


getMplusNames <- function(x, base = x, taken = character(0L)) {
  # Valid names are kept as is (unless they collide with names already taken),
  # whilst the rest are (shortened and) renamed. Mplus is case insensitive, so
  # uniqueness must be checked ignoring case.
  keep    <- !duplicated(x)
  base    <- base[keep] # NOTE: must be evaluated before `x` is modified
  x       <- x[keep]
  out     <- structure(character(length(x)), names = x)

  isValid <- x == base & isValidMplusName(x) &
    !duplicated(toupper(x)) & !toupper(x) %in% toupper(taken)

  out[isValid] <- x[isValid]
  taken <- c(taken, x[isValid])

  for (i in which(!isValid)) {
    clean <- stringr::str_remove_all(base[[i]], "[^A-Za-z0-9_]")
    if (!grepl("^[A-Za-z]", clean)) clean <- paste0("V", clean)

    candidate <- substr(clean, 1L, MPLUS_MAX_NAME_LENGTH)
    suffix    <- 0L
    while (!isValidMplusName(candidate) || toupper(candidate) %in% toupper(taken)) {
      suffix    <- suffix + 1L
      maxlen    <- MPLUS_MAX_NAME_LENGTH - nchar(suffix)
      candidate <- paste0(substr(clean, 1L, maxlen), suffix)
    }

    out[[i]] <- candidate
    taken    <- c(taken, candidate)
  }

  out
}


isValidMplusName <- function(x) {
  grepl("^[A-Za-z][A-Za-z0-9_]*$", x) &
    nchar(x) <= MPLUS_MAX_NAME_LENGTH &
    !toupper(x) %in% MPLUS_RESERVED_WORDS
}


specifyInteractionsMplus <- function(vars, names) {
  # Interaction terms are specified using `XWITH`, which requires at least one of
  # the variables to be latent. Interactions between observed variables are
  # computed as product terms using DEFINE instead. `XWITH` only takes two variables,
  # so higher order interactions are specified recursively (e.g., X:Z:W = (X:Z):W).
  XWITH    <- NULL
  DEFINE   <- NULL
  products <- NULL
  observed <- NULL
  intNames <- NULL

  getName <- function(xz) {
    allNames <- c(names, intNames)

    if (!xz %in% names(allNames)) {
      base     <- paste0(allNames[stringr::str_split_1(xz, ":")], collapse = "")
      intNames <<- c(intNames, getMplusNames(xz, base = base, taken = allNames))
    }

    c(names, intNames)[[xz]]
  }

  for (xz in vars$intTerms) {
    elems    <- stringr::str_split_1(xz, ":")
    isLatent <- elems %in% vars$latent
    obsElems <- elems[!isLatent]

    # Observed variables are multiplied into a single product term, which is then
    # (if there are any latent variables) interacted with the latent variables
    if (length(obsElems) > 1L) {
      product  <- paste0(obsElems, collapse = ":")
      products <- union(products, product)
      observed <- union(observed, obsElems)
      DEFINE   <- union(DEFINE, sprintf("%s = %s;", getName(product),
                                        paste0(names[obsElems], collapse = "*")))
      obsElems <- product
    }

    chain <- c(elems[isLatent], obsElems)

    if (length(chain) == 1L) {
      getName(xz)
      next
    }

    observed <- union(observed, setdiff(obsElems, products))

    for (j in seq_along(chain)[-1L]) {
      left  <- paste0(chain[seq_len(j - 1L)], collapse = ":")
      right <- chain[[j]]
      term  <- if (j == length(chain)) xz else paste0(chain[seq_len(j)], collapse = ":")

      XWITH <- union(XWITH, sprintf("%s | %s XWITH %s;",
                                    getName(term), getName(left), getName(right)))
    }
  }

  list(XWITH = XWITH, DEFINE = DEFINE, products = products,
       observed = observed, names = intNames)
}


getLevelVariablesMplus <- function(used, within = NULL, between = NULL) {
  if (!all(c("within", "between") %in% names(used)))
    return(list(within = within, between = between))

  if (is.null(within))  within  <- setdiff(used$within, used$between)
  if (is.null(between)) between <- setdiff(used$between, used$within)

  list(within = within, between = between)
}


mplusVariableList <- function(option, variables) {
  if (!length(variables))
    return(NULL)

  line <- sprintf("%s = %s;", option, paste0(variables, collapse = " "))
  paste0(strwrap(line, width = MPLUS_MAX_LINE_WIDTH, exdent = 2), collapse = "\n")
}


parTableToMplusModel <- function(parTable, names) {
  rename <- \(x) ifelse(x %in% names(names), names[x], x)

  lhs <- rename(parTable$lhs)
  rhs <- rename(parTable$rhs)
  op  <- parTable$op
  mod <- parTable$mod

  isRandom <- grepl("^rv\\(", mod)
  slopes   <- rename(stringr::str_match(mod, "^rv\\(\\s*[\"']?([^\"')]+)")[, 2L])

  modifier <- ifelse(mod == "", yes = "",
              ifelse(mod == "NA", yes = "*",
              ifelse(canBeNumeric(mod), yes = paste0("@", mod),
              ifelse(isRandom, yes = "",
                     no = paste0(" (", rename(mod), ")")))))

  statement <- character(NROW(parTable))
  for (i in seq_len(NROW(parTable))) {
    statement[[i]] <- switch(op[[i]],
      "=~" = sprintf("%s BY %s%s;", lhs[[i]], rhs[[i]], modifier[[i]]),
      "~"  = sprintf("%s ON %s%s;", lhs[[i]], rhs[[i]], modifier[[i]]),
      "~1" = sprintf("[%s%s];", lhs[[i]], modifier[[i]]),
      "|"  = sprintf("[%s$%s%s];", lhs[[i]], stringr::str_remove(rhs[[i]], "^t"), modifier[[i]]),
      "~~" = if (lhs[[i]] == rhs[[i]]) sprintf("%s%s;", lhs[[i]], modifier[[i]])
             else sprintf("%s WITH %s%s;", lhs[[i]], rhs[[i]], modifier[[i]]),
      mod_msg_stop(sprintf("Operator not supported for use in Mplus: `%s`", op[[i]]))
    )

    if (isRandom[[i]]) {
      statement[[i]] <- sprintf("%s | %s", slopes[[i]],
                                stringr::str_remove(statement[[i]], "\\*$"))
    }
  }

  statement
}


parTableToMplusConstraints <- function(constraints, names) {
  if (!NROW(constraints))
    return(NULL)

  labels <- names[names(names) != names]
  renameLabels <- function(expr) {
    for (label in names(labels)) {
      pattern <- sprintf("(?<![A-Za-z0-9_.])%s(?![A-Za-z0-9_.])",
                         stringr::str_escape(label))
      expr <- gsub(pattern, labels[[label]], expr, perl = TRUE)
    }

    stringr::str_replace_all(expr, stringr::fixed("^"), "**")
  }

  lhs <- renameLabels(constraints$lhs)
  rhs <- renameLabels(constraints$rhs)
  op  <- ifelse(constraints$op %in% c(":=", "=="), yes = "=", no = constraints$op)

  new <- lhs[constraints$op == ":="]
  NEW <- if (length(new)) sprintf("NEW(%s);", paste0(new, collapse = " ")) else NULL

  statements <- c(NEW, sprintf("%s %s %s;", lhs, op, rhs))
  paste0(strwrap(statements, width = MPLUS_MAX_LINE_WIDTH, exdent = 2), collapse = "\n")
}


specifyMplusAnalysis <- function(spec,
                                 estimator = "ml",
                                 type = NULL,
                                 algorithm = NULL,
                                 processors = 2,
                                 integration = 15,
                                 options = NULL) {
  isBayes <- grepl("bayes", estimator, ignore.case = TRUE)

  if (is.null(type)) {
    random <- spec$has.xwith || NROW(spec$slopes) > 0L
    base   <- if (spec$multilevel) "twolevel" else if (spec$has.cluster) "complex"

    type <- paste(c(base, if (random) "random"), collapse = " ")
    if (!nchar(type)) type <- "general"
  }

  # Numerical integration is not used by the Bayesian estimator
  if (is.null(algorithm) && spec$has.xwith && !isBayes)
    algorithm <- "integration"

  isIntegration <- !is.null(algorithm) && grepl("integration", algorithm, ignore.case = TRUE)

  defaults <- c(
    ESTIMATOR   = estimator,
    TYPE        = type,
    ALGORITHM   = algorithm,
    PROCESSORS  = processors,
    INTEGRATION = if (isIntegration) integration
  )

  # Additional (user specified) options, which take precedence over the defaults
  options <- unlist(options)
  mod_stopif(length(options) && (is.null(names(options)) || any(names(options) == "")),
             "`analysis.options` must be a named list!")
  if (length(options)) names(options) <- toupper(names(options))

  options <- c(defaults[!names(defaults) %in% names(options)], options)

  list(
    estimator = estimator,
    type      = type,
    algorithm = algorithm,
    ANALYSIS  = paste0(sprintf("%s = %s;", names(options), options), collapse = "\n")
  )
}


specifyMplusOutput <- function(spec, analysis, output.std = TRUE) {
  kway <- max(0L, lengths(stringr::str_split(spec$intTerms, ":")))

  # Mplus only supports standardized estimates for TYPE = TWOLEVEL RANDOM
  # when using a Bayesian estimator
  isTwoLevelRandom <- grepl("twolevel", analysis$type, ignore.case = TRUE) &&
    grepl("random", analysis$type, ignore.case = TRUE)
  isBayes <- grepl("bayes", analysis$estimator, ignore.case = TRUE)
  std     <- output.std && kway <= 2L && (!isTwoLevelRandom || isBayes)

  paste0(c("TECH1;", "TECH3;", if (std) "STANDARDIZED;"), collapse = "\n")
}


# ──────────────────────────────────────────────────────────────────────────────
# Running Mplus
# ──────────────────────────────────────────────────────────────────────────────
runMplusModel <- function(input, cleanup = FALSE) {
  fprefix <- getMplusFilePrefix()
  if (cleanup)
    on.exit(cleanupMplusFiles(fprefix), add = TRUE)

  inp <- paste0(fprefix, ".inp")
  out <- paste0(fprefix, ".out")

  # Write the input and data files, and run the model separately, such that we
  # only read the parts of the output we need
  MplusAutomation::mplusModeler(input, modelout = inp, run = 0L,
                                writeData = "always", hashfilename = FALSE)
  MplusAutomation::runModels(inp, logFile = NULL, quiet = TRUE)

  mod_stopif(!file.exists(out), "Unable to find `Mplus` output file!")

  results <- MplusAutomation::readModels(
    out, quiet = TRUE,
    what = c("input", "warn_err", "summaries", "parameters", "tech1", "tech3")
  )

  # Some errors (e.g., a non-positive definite information matrix) still
  # produce parameter estimates, in which case we only warn
  errors <- vapply(results$errors, FUN.VALUE = character(1L),
                   FUN = \(e) paste0(e, collapse = " "))
  failed <- is.null(results$parameters$unstandardized)

  if (failed && length(errors)) {
    mod_msg_stop("`Mplus` returned the following error(s):", errors)
  } else if (failed) {
    mod_msg_stop("Unable to retrieve parameter estimates from `Mplus`!",
                 sprintf("Check the output file (%s) for errors.", out))
  } else if (length(errors)) {
    mod_msg_warn("`Mplus` returned the following error(s):", errors)
  }

  input$results <- results
  input
}


mplusFilePrefixExists <- function(fprefix) {
  any(file.exists(paste0(fprefix, c(".inp", ".dat", ".out"))))
}


cleanupMplusFiles <- function(fprefix) {
  files <- paste0(fprefix, c(".inp", ".dat", ".out"))
  unlink(files[file.exists(files)])
}


getMplusFilePrefix <- function(max.iter = 20, base = "mpresults") {
  fprefix <- paste0(base, generateRandomCharId(n = 12))

  for (i in seq_len(max.iter)) {
    if (!mplusFilePrefixExists(fprefix))
      return(fprefix)

    fprefix <- paste0(fprefix, generateRandomCharId(2))
  }

  mod_msg_warn(
    "Unable to create a unique name for Mplus files!",
    "Previous results might get overwritten..."
  )

  fprefix
}


generateRandomCharId <- function(n = 36) {
  # use only lower letters, as Mplus sometimes seem to convert upper case
  # characters in the .inp file to lower characters in the .out file. In
  # particular this seems to happen on Mplus 8.11 on Windows
  chars <- c(letters, as.character(0:9))
  paste0(sample(chars, size = n, replace = TRUE), collapse = "")
}


# ──────────────────────────────────────────────────────────────────────────────
# Reading Mplus Results
# ──────────────────────────────────────────────────────────────────────────────
mplusParTable <- function(params, spec) {
  # Convert a parameter table from MplusAutomation (e.g.,
  # `results$parameters$unstandardized`) into a modsem/lavaan parameter table
  if (is.null(params) || !NROW(params))
    return(NULL)

  header <- params$paramHeader
  param  <- params$param

  isBY   <- endsWith(header, ".BY")
  isON   <- endsWith(header, ".ON")
  isWITH <- endsWith(header, ".WITH")
  isMean <- header %in% c("Means", "Intercepts")
  isVar  <- header %in% c("Variances", "Residual.Variances")
  isThr  <- header == "Thresholds"
  isNew  <- header == "New.Additional.Parameters"

  known <- isBY | isON | isWITH | isMean | isVar | isThr | isNew
  mod_warnif(!all(known),
             sprintf("Ignoring unrecognized parameter(s) from Mplus: %s",
                     paste0(unique(header[!known]), collapse = ", ")))

  # Standardized random slopes are reported as e.g., `S|FW ON X1`, and are
  # mapped to the (within level) effect they belong to (i.e., `fw ~ x1`)
  lhs <- stringr::str_remove(header, "\\.(BY|ON|WITH)$")
  lhs <- stringr::str_remove(lhs, "^.*\\|")
  lhs[isMean | isVar | isNew] <- param[isMean | isVar | isNew]
  lhs[isThr] <- stringr::str_split_i(param[isThr], "\\$", i = 1L)

  rhs <- param
  rhs[isMean | isNew] <- ""
  rhs[isThr] <- paste0("t", stringr::str_split_i(param[isThr], "\\$", i = 2L))

  op <- rep(NA_character_, length(header))
  op[isBY]   <- "=~"
  op[isON]   <- "~"
  op[isWITH | isVar] <- "~~"
  op[isMean] <- "~1"
  op[isThr]  <- "|"
  op[isNew]  <- ":="

  level <- if ("BetweenWithin" %in% names(params)) params$BetweenWithin else rep("", length(op))
  level <- tolower(level)

  na999 <- \(x) ifelse(x == 999, yes = NA, no = x)

  if ("posterior_sd" %in% names(params)) {
    # Bayesian estimates: posterior standard deviations, one-tailed p-values,
    # and credible intervals. Fixed parameters have a posterior sd of 0.
    fixed    <- params$posterior_sd == 0
    se       <- ifelse(fixed, yes = NA, no = params$posterior_sd)
    z        <- rep(NA_real_, length(se))
    ci.lower <- ifelse(fixed, yes = NA, no = params$lower_2.5ci)
    ci.upper <- ifelse(fixed, yes = NA, no = params$upper_2.5ci)

  } else {
    # Mplus reports fixed parameters with a standard error of 0, and a p-value of 999
    fixed    <- params$se == 0 & params$pval == 999
    se       <- ifelse(fixed, yes = NA, no = params$se)
    z        <- ifelse(fixed, yes = NA, no = na999(params$est_se))
    ci.lower <- params$est - CI_WIDTH * se
    ci.upper <- params$est + CI_WIDTH * se
  }

  parTable <- data.frame(
    lhs       = lhs,
    op        = op,
    rhs       = rhs,
    level     = level,
    label     = "",
    est       = params$est,
    std.error = se,
    z.value   = z,
    p.value   = ifelse(fixed, yes = NA, no = na999(params$pval)),
    ci.lower  = ci.lower,
    ci.upper  = ci.upper,
    stringsAsFactors = FALSE
  )[known, , drop = FALSE]

  parTable <- renameParTableFromMplus(parTable, spec = spec)
  rownames(parTable) <- NULL

  parTable
}


renameParTableFromMplus <- function(parTable, spec) {
  # Mplus converts all names to upper case
  names <- structure(names(spec$names), names = toupper(spec$names))
  rename <- \(x) ifelse(toupper(x) %in% names(names), names[toupper(x)], x)

  isNew <- parTable$op == ":="
  parTable$lhs <- rename(parTable$lhs)
  parTable$rhs[!isNew] <- rename(parTable$rhs[!isNew])

  # Defined parameters
  custom <- spec$constraints[spec$constraints$op == ":=", , drop = FALSE]
  parTable$rhs[isNew]   <- trimws(custom$rhs[match(parTable$lhs[isNew], custom$lhs)])
  parTable$label[isNew] <- parTable$lhs[isNew]

  # Labels
  input   <- spec$parTable
  isLabel <- input$mod != "" & !canBeNumeric(input$mod) &
    !grepl("^rv\\(", input$mod) & input$mod != "NA"
  input   <- input[isLabel, , drop = FALSE]

  key <- \(pt) paramKeysMplus(pt$lhs, pt$op, pt$rhs, pt$level)
  match <- match(key(parTable), key(input))
  parTable$label[!is.na(match)] <- input$mod[match[!is.na(match)]]

  # Levels (defined parameters are not part of any level)
  if (spec$multilevel) {
    parTable$level <- ifelse(parTable$level == "between", yes = 2L, no = 1L)
    parTable$level[isNew] <- 0L
  } else parTable$level <- NULL

  parTable
}


paramKeysMplus <- function(lhs, op, rhs, level) {
  n   <- max(length(lhs), length(op), length(rhs), length(level))
  lhs <- rep_len(lhs, n)
  rhs <- rep_len(rhs, n)

  # Covariances are symmetric, so the order of lhs and rhs should not matter
  swap <- op == "~~" & lhs > rhs
  tmp  <- lhs[swap]
  lhs[swap] <- rhs[swap]
  rhs[swap] <- tmp

  paste(toupper(lhs), op, toupper(rhs), tolower(level), sep = "\r")
}


addMplusParamIds <- function(parTable, TECH1, spec) {
  # Get the parameter numbers in TECH1 (which define the order of parameters
  # in the vcov matrix in TECH3) for each parameter in the parameter table
  ids <- getMplusParamIdsTECH1(TECH1)

  if (!NROW(ids)) {
    mod_msg_warn("Unable to retrieve parameter specification (TECH1) from `Mplus`!")
    parTable$id <- NA_integer_
    return(parTable)
  }

  # Defined parameters are not part of any level, and have no rhs in TECH1
  isNew <- parTable$op == ":="
  level <- if (spec$multilevel) c("", "within", "between")[parTable$level + 1L] else ""
  rhs   <- ifelse(isNew, yes = "", no = renameParTableToMplus(parTable$rhs, spec))
  keys  <- paramKeysMplus(renameParTableToMplus(parTable$lhs, spec), parTable$op,
                          rhs, level)

  # A parameter in a TECH1 matrix can correspond to different parameters in the
  # parameter table. E.g., an element in LAMBDA can be both `f =~ y` and `y ~ f`,
  # whilst an element in BETA can be both `f ~ x` and `x =~ f` (higher order models)
  parTable$id <- NA_integer_
  for (alt in c("", ".alt")) {
    altKeys <- paramKeysMplus(ids[[paste0("lhs", alt)]], ids[[paste0("op", alt)]],
                              ids[[paste0("rhs", alt)]], ids$level)
    match   <- match(keys, altKeys)
    fill    <- is.na(parTable$id) & !is.na(match)
    parTable$id[fill] <- ids$id[match[fill]]
  }

  missing <- setdiff(ids$id, parTable$id)
  mod_warnif(length(missing),
             sprintf("Unable to match parameter(s) %s in TECH1 to parameter table!",
                     paste0(missing, collapse = ", ")))

  parTable
}


renameParTableToMplus <- function(x, spec) {
  ifelse(x %in% names(spec$names), spec$names[x], x)
}


getMplusParamIdsTECH1 <- function(TECH1) {
  spec <- TECH1$parameterSpecification

  if (is.null(spec) || !length(spec))
    return(NULL)

  custom <- spec[["THE.ADDITIONAL.PARAMETERS"]]$new_additional
  spec   <- spec[names(spec) != "THE.ADDITIONAL.PARAMETERS"]

  blocks <- if (all(c("WITHIN", "BETWEEN") %in% names(spec))) {
    spec[c("WITHIN", "BETWEEN")]
  } else if (!is.null(spec$lambda) || !is.null(spec$nu) || !is.null(spec$beta)) {
    list(spec)
  } else spec[1L]

  levels <- if (length(blocks) > 1L) c("within", "between") else ""

  ids <- lapply(seq_along(blocks), FUN = \(i) {
    getMplusParamIdsBlock(blocks[[i]], level = levels[[i]])
  })

  rbind(do.call(rbind, ids), getMplusParamIdsMatrix(
    custom, op = ":=", lhs = \(row, col) col, rhs = \(row, col) "",
    level = ""
  ))
}


getMplusParamIdsBlock <- function(block, level) {
  vec   <- \(row, col) col
  empty <- \(row, col) ""
  row   <- \(row, col) row
  thrLhs <- \(row, col) stringr::str_split_i(col, "\\$", i = 1L)
  thrRhs <- \(row, col) paste0("t", stringr::str_split_i(col, "\\$", i = 2L))

  rbind(
    getMplusParamIdsMatrix(block$tau,   op = "|",  lhs = thrLhs, rhs = thrRhs, level = level),
    getMplusParamIdsMatrix(block$nu,    op = "~1", lhs = vec, rhs = empty, level = level),
    getMplusParamIdsMatrix(block$alpha, op = "~1", lhs = vec, rhs = empty, level = level),
    getMplusParamIdsMatrix(block$lambda, op = "=~", lhs = vec, rhs = row, level = level,
                           op.alt = "~", lhs.alt = row, rhs.alt = vec),
    getMplusParamIdsMatrix(block$beta,  op = "~", lhs = row, rhs = vec, level = level,
                           op.alt = "=~", lhs.alt = vec, rhs.alt = row),
    getMplusParamIdsMatrix(block$gamma, op = "~", lhs = row, rhs = vec, level = level),
    getMplusParamIdsMatrix(block$psi,   op = "~~", lhs = row, rhs = vec, level = level),
    getMplusParamIdsMatrix(block$theta, op = "~~", lhs = row, rhs = vec, level = level)
  )
}


getMplusParamIdsMatrix <- function(M, op, lhs, rhs, level,
                                   op.alt = op, lhs.alt = lhs, rhs.alt = rhs) {
  if (is.null(M) || !length(M))
    return(NULL)

  M     <- as.matrix(M)
  index <- which(!is.na(M) & M > 0, arr.ind = TRUE)

  if (!NROW(index))
    return(NULL)

  rows <- rownames(M)[index[, 1L]]
  cols <- colnames(M)[index[, 2L]]
  if (is.null(rows)) rows <- rep("", NROW(index))

  data.frame(
    id      = as.integer(M[index]),
    lhs     = lhs(rows, cols),
    op      = op,
    rhs     = rhs(rows, cols),
    lhs.alt = lhs.alt(rows, cols),
    op.alt  = op.alt,
    rhs.alt = rhs.alt(rows, cols),
    level   = level,
    stringsAsFactors = FALSE
  )
}


getFreeParamsMplus <- function(parTable) {
  # The (unique) free parameters in the parameter table, with their row in the
  # parameter table, and their parameter number (id) in Mplus
  free <- which(!is.na(parTable$id) & !duplicated(parTable$id))
  free <- free[order(parTable$id[free])] # Sort by the order used by Mplus

  data.frame(
    row   = free,
    id    = parTable$id[free],
    label = getParTableLabels(parTable[free, , drop = FALSE], labelCol = "label")
  )
}


removeParamIdsMplus <- function(parTable) {
  parTable[names(parTable) != "id"]
}


getVcovMplus <- function(model, free) {
  labels <- free$label

  vcov <- tryCatch({
    V <- model$results$tech3$paramCov
    V[upper.tri(V)] <- t(V)[upper.tri(V)]
    V[free$id, free$id, drop = FALSE]
  }, error = function(e) {
    mod_msg_warn(paste0("Unable to retrive `tech3` from `Mplus` results\n",
                        "Message: ", conditionMessage(e)))
    matrix(NA_real_, nrow = NROW(free), ncol = NROW(free))
  })

  dimnames(vcov) <- list(labels, labels)
  vcov
}
