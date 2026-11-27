#' Deprecated functions
#'
#' @description
#' These functions are deprecated:
#'
#' * `havefftw()` is replaced by [has_fftw()].
#' * `havealglib()` always returns `FALSE`, since chebpol no longer supports ALGLIB.
#'
#' @keywords internal
#' @name chebpol-deprecated
NULL

#' @rdname chebpol-deprecated
#' @export
havefftw <- function() {
  .Deprecated("has_fftw", package = "chebpol")
  has_fftw()
}

#' @rdname chebpol-deprecated
#' @export
havealglib <- function() {
  .Deprecated(msg = "'havealglib' is deprecated. chebpol no longer supports ALGLIB.")
  FALSE
}
