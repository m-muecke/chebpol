.onAttach <- function(libname, pkgname) {
  if (!havefftw()) {
    packageStartupMessage(
      "*** ",
      pkgname,
      ": FFTW not used.\n*** You should install it from https://fftw.org\n*** or check if your OS-distribution provides it, and recompile."
    )
  }
}

.onLoad <- function(libname, pkgname) {
  thr <- as.integer(Sys.getenv("CHEBPOL_THREADS"))
  if (is.na(thr)) {
    thr <- 1L
  }
  options(chebpol.threads = thr)
}
