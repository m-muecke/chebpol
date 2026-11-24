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
  threads <- as.integer(Sys.getenv("CHEBPOL_THREADS", unset = "1"))
  options(chebpol.threads = getOption("chebpol.threads", if (is.na(threads)) 1L else threads))
}
