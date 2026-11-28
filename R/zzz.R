.onLoad <- function(libname, pkgname) {
  threads <- as.integer(Sys.getenv("CHEBPOL_THREADS", unset = "1"))
  options(chebpol.threads = getOption("chebpol.threads", if (is.na(threads)) 1L else threads))
}
