getMethodSettingsDA <- function(method, args = NULL) {

    settings <- list(
        lms = list(verbose = interactive(),
                   optimize = TRUE,
                   # Resolved below by `quadNodeDefault()` from `integration`,
                   # `adaptive` and `lms.backend`. This is the value for the
                   # default legacy `gh` + `quasi` combination.
                   nodes = 24,
                   convergence.abs = 1e-4,
                   convergence.rel = 1e-10,
                   optimizer = "nlminb",
                   center.data = FALSE,
                   standardize.data = FALSE,
                   standardize.out = FALSE,
                   standardize = FALSE,
                   mean.observed = TRUE,
                   double = FALSE,
                   calc.se = TRUE,
                   FIM = "observed",
                   OFIM.hessian = TRUE,
                   EFIM.S = 100,
                   EFIM.parametric = TRUE,
                   robust.se = FALSE,
                   R.max = 1e5,
                   max.iter = 500,
                   max.step = 1,
                   epsilon = 1e-6,
                   quad.range = Inf,
                   rect.range = 5,
                   lms.backend = "legacy",
                   integration = "gh",
                   adaptive = "quasi",
                   adaptive.quad.tol = 1e-12,
                   adaptive.frequency = 3,
                   n.threads = NULL,
                   algorithm = "EMA",
                   em.control = list(),
                   missing = "listwise",
                   orthogonal.x = FALSE,
                   orthogonal.y = FALSE,
                   auto.fix.first = TRUE,
                   auto.fix.single = TRUE,
                   auto.split.syntax = FALSE,
                   cr1s = FALSE,
                   compress.data = FALSE,
                   group = NULL,
                   sampling.weights = NULL,
                   sampling.weights.normalization = "total",
                   fix.composite.var = TRUE
        ),
        qml = list(verbose = interactive(),
                   optimize = TRUE,
                   nodes = 0,
                   convergence.rel = 1e-6,
                   convergence.abs = NULL, # not relevant
                   optimizer = "nlminb",
                   center.data = FALSE,
                   standardize = FALSE,
                   standardize.data = FALSE,
                   standardize.out = FALSE,
                   mean.observed = TRUE,
                   double = FALSE,
                   calc.se = TRUE,
                   FIM = "observed",
                   OFIM.hessian = TRUE,
                   EFIM.S = 100,
                   EFIM.parametric = TRUE,
                   robust.se = FALSE,
                   R.max = 1e5,
                   max.iter = 500,
                   max.step = NULL,
                   epsilon = 1e-8,
                   quad.range = Inf,
                   rect.range = 5,
                   lms.backend = "legacy",
                   integration = "gh",
                   adaptive = "none",
                   adaptive.quad.tol = NULL,
                   n.threads = NULL,
                   adaptive.frequency = NULL,
                   em.control = NULL,
                   algorithm = NULL,
                   missing = "listwise",
                   orthogonal.x = FALSE,
                   orthogonal.y = FALSE,
                   auto.fix.first = TRUE,
                   auto.fix.single = TRUE,
                   auto.split.syntax = TRUE,
                   cr1s = FALSE,
                   compress.data = FALSE,
                   group = NULL,
                   sampling.weights = NULL,
                   sampling.weights.normalization = "total",
                   fix.composite.var = FALSE
        )
    )

    if (is.null(args)) return(settings[method])

    settingNames <- unique(unlist(lapply(settings, FUN = names)))
    args <- args[settingNames]
    isMissing <- vapply(args, FUN.VALUE = logical(1L), FUN = is.null)
    missingArgs <- settingNames[isMissing]

    mod_stopif(!method %in% names(settings), "Unrecognized method")

    args.out <- c(settings[[method]][missingArgs], args[!isMissing])

    if (method == "lms") {
      args.out$integration <- match.arg(args.out$integration, INTEGRATION_TYPES)
      args.out$adaptive    <- match.arg(args.out$adaptive, ADAPTIVE_TYPES)

      # How many nodes are sensible depends on both axes: a per-observation
      # rule needs far fewer points than a shared one, and Monte-Carlo counts
      # total draws rather than nodes per dimension.
      if (is.null(args$nodes))
        args.out$nodes <- quadNodeDefault(args.out$integration,
                                          args.out$adaptive,
                                          args.out$lms.backend %||% "legacy")

      # The graph backend costs one N x Q traversal per objective evaluation,
      # so its runtime is linear in the number of rows. Collapsing duplicate
      # rows is exact (see `compressData()`), and on the categorical data this
      # backend targets it removes a large fraction of them -- which is also
      # what Mplus does, so it keeps the two comparable. The legacy LMS and
      # QML backends are left alone: they are not row-bound in the same way,
      # and their weighted paths are less well exercised.
      if (is.null(args$compress.data))
        args.out$compress.data <-
          identical(args.out$lms.backend %||% "legacy", "graph")
    }

    args.out$standardize.data <-
      args.out$standardize || args.out$standardize.data
    args.out$standardize.out <-
      args.out$standardize || args.out$standardize.out
    args.out$mean.observed <-
      !args.out$standardize && args.out$mean.observed
    args.out$OFIM.hessian <-
      args.out$OFIM.hessian && !args.out$robust.se
    args.out$center.data <- !args.out$standardize.data &&
      !args.out$mean.observed

    if (is.null(args.out$group) &&
        tolower(args.out$sampling.weights.normalization) == "group")
      args.out$sampling.weights.normalization <- "total"

    args.out
}



#' default arguments fro LMS and QML approach
#'
#' @param method which method to get the settings for
#' @return list
#' @export
#' @description
#' This function returns the default settings for the LMS and QML approach.
#' @examples
#' library(modsem)
#' default_settings_da()
default_settings_da <- function(method = c("lms", "qml")) {
  getMethodSettingsDA(method = method, args = NULL)
}
