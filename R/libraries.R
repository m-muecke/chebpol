#' Check which optional libraries chebpol was built with
#'
#' chebpol can be compiled with or without three optional C libraries, each of which enables part of the package:
#'
#' * FFTW computes Chebyshev coefficients with a fast DCT. Without it, `chebcoef()` and `ipol()` with
#'   `method = "chebyshev"` use a matrix product, which is not feasible for high-degree polynomials: a 100x100x100
#'   approximation is possible, but not a one-dimensional 1000000.
#' * ALGLIB provides compact support radial basis function interpolation of scattered data, `ipol()` with
#'   `method = "crbf"`.
#' * GSL solves for the degree of the stalker splines on a non-uniform grid. Without it, `ipol()` with
#'   `method = "stalker"` uses degree 2 there and warns, and `method = "oldstalker"` with `k = NA` stops with an error.
#'
#' When chebpol is installed from source on Linux or macOS, FFTW and GSL are used if pkg-config finds them. To build
#' without one, set the environment variable `CHEBPOL_DISABLE_FFTW` or `CHEBPOL_DISABLE_GSL` to a non-empty value
#' before installing, or pass `--disable-fftw` or `--disable-gsl` in `configure.args`.
#'
#' @return `TRUE` if chebpol was built with the library, otherwise `FALSE`.
#' @export
#' @examples
#' has_fftw()
#' has_alglib()
#' has_gsl()
has_fftw <- function() .Call(c_has_fftw)

#' @rdname has_fftw
#' @export
has_alglib <- function() .Call(c_has_alglib)

#' @rdname has_fftw
#' @export
has_gsl <- function() .Call(c_has_gsl)
