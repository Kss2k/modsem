#ifndef THREAD_SETTER_H
#define THREAD_SETTER_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Rcpp.h>

// RAII wrapper: sets OMP thread count on construction, restores on destruction.
struct ThreadSetter {
  int old_threads = 1;

  explicit ThreadSetter(const int ncores) {
#ifdef _OPENMP
    if (ncores <= 0) Rcpp::stop("ncores must be positive");
    old_threads = omp_get_max_threads();
    omp_set_num_threads(ncores);
#endif
  }

  ~ThreadSetter() {
#ifdef _OPENMP
    omp_set_num_threads(old_threads);
#endif
  }
};

#endif // THREAD_SETTER_H
