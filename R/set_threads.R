ThreadEnv <- rlang::env(n.threads  = NULL,
                        n.blas     = NULL,
                        n.blas.old = NULL)


blas_set_num_threads <- function(n, .default = 1L) {
  tryCatch(
    purrr::quietly(RhpcBLASctl::blas_set_num_threads)(n),
    error = function(e) {
      warning2("Could not set threads for OpenBLAS!\nMessage: ",
               conditionMessage(e))

      .default # return default
    }
  )
}


blas_get_num_procs <- function(.default = 1L) {
  tryCatch(
    purrr::quietly(RhpcBLASctl::blas_get_num_procs)(),
    error = function(e) {
      warning2("Could not get processors for OpenBLAS!\nMessage: ",
               conditionMessage(e))

      .default # return default
    }
  )
}


setThreadEnv <- function(n.threads, n.blas = NULL) {
  ThreadEnv$n.threads <- n.threads

  if (!is.null(n.blas)) {
    ThreadEnv$n.blas.old <- blas_get_num_procs()
    ThreadEnv$n.blas     <- n.blas

    blas_set_num_threads(n.blas)

  } else if (!is.null(ThreadEnv$n.blas.old)) {
    # If we're not setting threads for BLAS, we should check if we need to
    # reset them back to some previous value

    blas_set_num_threads(ThreadEnv$n.blas.old)

    ThreadEnv$n.blas     <- NULL
    ThreadEnv$n.blas.old <- NULL
  }
}


resetThreadEnv <- function() {
  setThreadEnv(n.threads = NULL, n.blas = NULL)
}


setThreads <- function(n, n.blas = 1L) {
  resetThreadEnv() # reset

  if (is.null(n))
    k <- getDefaultThreads()

  else if (is.numeric(n))
    k <- getThreadsN(n)
  else if (is.character(n)) {
    k <- switch(n,
                default = getDefaultThreads(),
                max     = getMaxThreads(),
                min     = getMinThreads(),
                stop2("Invalid string specifying number of threads"))

  } else stop2("Invalid number of threads, must be integer, NULL, or character")

  setThreadEnv(n.threads = k, n.blas = n.blas)
}


resetThreads <- function() {
  resetThreadEnv()
}


getDefaultThreads <- function() {
  defaultCRAN <- 2
  mc.cores    <- options("mc.cores")[[1L]]
  default     <- if (!is.null(mc.cores)) mc.cores else defaultCRAN

  getThreadsN(default)
}


getMaxThreads <- function() {
  parallel::detectCores()
}


getMinThreads <- function() {
  getThreadsN(n = 1)
}


getThreadsN <- function(n) {
  ncores <- parallel::detectCores()
  ifelse(n >= ncores, ncores, n)
}
