.onAttach <- function(libname, pkgname) {
  if (!havefftw()) {
    packageStartupMessage(
      "*** ",
      pkgname,
      ": FFTW not used.\n*** You should install it from http://fftw.org\n*** or check if your OS-distribution provides it, and recompile."
    )
  }
}

.onLoad <- function(libname, pkgname) {
  if (is.na(thr <- as.integer(Sys.getenv("CHEBPOL_THREADS")))) {
    thr <- 1L
  }
  options(chebpol.threads = thr)
}

deprecated <- function(fun, ...) {
  realfun <- call(":::", quote(chebpol), as.name(paste0(fun, ".real")))
  .Deprecated("ipol", "chebpol", old = fun)
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- realfun
  mc <- mc[-2]
  eval.parent(mc)
}
