#' Deprecated functions
#'
#' @description
#' These functions are deprecated:
#'
#' * `havefftw()` is replaced by [has_fftw()].
#' * `havealglib()` is replaced by [has_alglib()].
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
  .Deprecated("has_alglib", package = "chebpol")
  has_alglib()
}
