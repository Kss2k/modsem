devtools::load_all()

m1 <- "
# Outer Model
  X =~ x1
  X =~ x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

# Inner model
  Y ~ X + Z + X:Z
"

est1 <- modsem(m1, data = oneInt)
simple_slopes(x = "X", z = "Z", y = "Y", model = est1)
plot_interaction(x = "X", z = "Z", y = "Y",
                 vals_z = c(1, 0), model = est1, greyscale = TRUE)
plot_interaction(x = "X", z = "Z", y = "Y", xz = "X:Z",
                 vals_z = c(1, 0), model = est1, ci_type = "prediction")
plot_jn(x = "X", z = "Z", y = "Y", model = est1, greyscale = TRUE)

plot_surface(x = "X", z = "Z", y = "Y", model = est1, colorscale = "Greys", grid = TRUE, grid_color = "black")
# check input length validation
checkInputValidation <- function(x.default = "X",
                                 z.default = "Z",
                                 y.default = "Y",
                                 xz.default = NULL,
                                 FUN = simple_slopes,
                                 debug = FALSE,
                                 ...) {
  params <- list(x  = x.default,
                 z  = z.default,
                 y  = y.default,
                 xz = xz.default)

  for (par in names(params)) {
    if (is.null(params[[par]]))
      next

    if (debug) {
      printf("  |> Checking %s, xz.default = %s\n", par,
             ifelse(is.null(xz.default), yes = "NULL", no = xz.default))
    }

    argsLong <- c(params, list(...))
    argsLong[[par]] <- c(argsLong[[par]], argsLong[[par]])

    testthat::expect_error(
      do.call(FUN, argsLong),
      regexp = sprintf("%s must be of length 1.*", par)
    )

    if (!par != "xz") { # xz is not validated the same way...
      argsWrong <- c(params, list(...))
      argsWrong[[par]] <- "I_DO_NOT_EXIST"
      testthat::expect_error(
         do.call(FUN, argsWrong),
         regexp = sprintf("Unrecognized variable: I_DO_NOT_EXIST\\!")
      )
    }
  }
}


debug <- TRUE
funcs <- list(simple_slopes, plot_interaction, plot_jn, plot_surface)
args  <- list(list(model = est1, vals_z = c(0, 1)),
              list(model = est1, vals_z = c(0, 1)),
              list(model = est1), list(model = est1))
for (i in seq_along(funcs)) {
  FUN <- funcs[[i]]
  if (debug) printf("Checking func %d...\n", i)
  for (xz_null in c(TRUE, FALSE)) {
    args_i <- c(args[[i]], list(FUN = FUN, debug = debug,
                                xz.default = if (xz_null) NULL else "X:Z"))

    do.call(checkInputValidation, args_i)
  }
}
