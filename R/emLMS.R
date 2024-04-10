
#' Em Lms
#'
#' @param model model object
#' @param verbose
#' @param convergence
#' @param maxiter
#' @param maxstep
#' @param max.singleClass
#' @param negHessian
#' @param breakOnLogLikIncrease
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
emLms <- function(model, verbose = FALSE,
               convergence = 1e-02, maxiter = 500,
               maxstep = 1, negHessian = TRUE,
               breakOnLogLikIncrease = FALSE,
               ...) {
  data <- model$data
  if (anyNA(data)) stop("Remove or replace missing values from data")

  logLikRet <- NULL
  logLikNew <- 0
  logLikOld <- Inf
  iterations <- 0     
  thetaNew <- model$theta
  nLogIncreased <- 0
  precision <- FALSE
  run <- TRUE
  while(run) { # as long as no convergence is reached
    if (logLikNew - logLikOld > 0) {
      if (breakOnLogLikIncrease) {
        warning("Loglikelihood is increasing. EM algorithm will be stopped.")
        logLikNew <- logLikOld
        thetaNew <- thetaOld
        break
      }
      nLogIncreased <- nLogIncreased + 1
      if (nLogIncreased > 5) {
        precision <- TRUE 
        maxstep <- 1
      }
    }


    # Update loglikelihood
    logLikOld <- logLikNew
    thetaOld <- thetaNew

    # E-step
    P <- estepLms(model = model, theta = thetaOld, dat = data, 
                  precision = precision, ...)
    # M-step
    mstep <- mstepLms(model = model, P = P, dat = data,
                      theta = thetaOld, maxstep = maxstep, ...,
                      precision = precision)

    logLikNew   <- mstep$objective
    logLikRet   <- c(logLikRet, logLikNew)
    thetaNew  <- unlist(mstep$par)
    iterations <- iterations + 1

    if (verbose) {
      cat(sprintf("EM: Iteration = %5d, LogLik = %11.2f, Change = %10.3f\n",
            iterations, -logLikNew, logLikOld - logLikNew))
    }
    if(iterations >= maxiter){
      warning("Maximum number of iterations was reached. ",
              "EM algorithm might not have converged.")
      break
    }
    if (abs(logLikOld - logLikNew) < convergence) run <- FALSE
  }
  final <- mstepLms(model = model, P = P, dat = data,
                    theta = thetaNew, negHessian = negHessian,
                    maxstep = maxstep, ...)
  coefficients <- final$par
  finalModel <- fillModel(model, coefficients, fillPhi = TRUE)


  # convergence of em
  if (iterations == maxiter) em_convergence <- FALSE else em_convergence <- TRUE


  info <- model$info
  info$iterations <- iterations

  out <- list(model = finalModel, coefficients=coefficients,
              objective=-final$objective, em.convergence=em_convergence,
              negHessian=final$hessian, loglikelihoods=-logLikRet, info=info)

  class(out) <- "modsemLMS"
  out
}
