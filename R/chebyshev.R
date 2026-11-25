# Chebyshev transformation.  I.e. coefficients for given function values in the knots.

# The Chebyshev knots of order n on an interval
chebknots1 <- function(n, interval = NULL) {
  kn <- cospi((seq_len(n) - 0.5) / n)
  if (is.null(interval)) {
    return(kn)
  }
  kn * diff(interval) / 2 + mean(interval)
}

#' Create a Chebyshev-grid
#'
#' Create a Chebyshev grid on a hypercube.
#'
#' If `intervals` is not provided, it is assumed that the domain of the
#' function in each dimension is `[-1, 1]`.  Thus, standard Chebyshev knots are
#' produced.  If `dims` is of length 1, `intervals` may be a vector
#' of length 2 rather than a list with a vector of length 2.
#'
#' @param dims The number of grid-points in each dimension.  For
#' Chebyshev-polynomial of degree `dims-1`.
#' @param intervals A list of vectors of length 2.  The lower and upper bounds
#' of the hypercube.
#' @return An array of dimension `dims`.  The Chebyshev grid-points.
#' @export
#' @examples
#' # Standard knots for degree 3
#' chebknots(4)
#' # Knots in the interval [2,3] for degree 3
#' chebknots(4, intervals = c(2, 3))
#' # Multivariate knots
#' chebknots(c(x = 3, y = 4, z = 3))
#' # Multivariate grid
#' head(expand.grid(chebknots(c(x = 3, y = 4, z = 5), list(c(1, 3), c(4, 6), c(800, 900)))))
chebknots <- function(dims, intervals = NULL) {
  if (is.null(intervals)) {
    res <- lapply(dims, chebknots1)
  } else {
    if (length(dims) == 1L && is.numeric(intervals)) {
      intervals <- list(intervals)
    }
    if (!is.list(intervals)) {
      stop("`intervals` must be a list with one length-2 vector per dimension.")
    }
    res <- Map(chebknots1, dims, intervals)
  }

  res
}

#' Compute Chebyshev-coefficients given values on a Chebyshev grid
#'
#' Compute the multivariate Chebyshev-coefficients, given values on a Chebyshev
#' grid.
#'
#' If `val` has no `dim`-attribute, it is assumed to be
#' one-dimensional of length the length of `val`.
#'
#' If chebpol was compiled without FFTW, running
#' `chebcoef` on large grids may be slow and memory-demanding.
#'
#' @param val An `array` of function values on a Chebyshev grid.  The
#' `dim`-attribute must be appropriately set.  If not set, it is assumed
#' to be one-dimensional.
#' @param dct Logical. Since the Chebyshev coefficients are closely related to
#' the DCT-II transform of `val`, the non-normalized real-even DCT-II
#' coefficients may be retrieved instead. I.e. those from FFTW_REDFT10 in each
#' dimension.  This is not used anywhere in the package, it is merely provided
#' as a convenience for those who might need it.
#' @return An array of Chebyshev-coefficients for an interpolating
#' Chebyshev-polynomial.
#' @seealso [havefftw()]
#' @export
#' @examples
#' # Coefficients for a 2x3x4 grid
#' a <- array(rnorm(24), dim = c(2, 3, 4))
#' chebcoef(a)
chebcoef <- function(val, dct = FALSE) {
  val <- as.array(val)
  if (is.integer(val) || is.logical(val)) {
    storage.mode(val) <- "double"
  }
  .Call(C_chebcoef, val, dct, getOption("chebpol.threads"))
}

#' Evaluate a Chebyshev interpolation in a point
#'
#' Given Chebyshev coefficients, evaluate the interpolation in a point.
#'
#' @param x The point to evaluate.
#' @param coef The Chebyshev coefficients. Typically from a call to
#' [chebcoef()], possibly modified. A vector without a `dim` attribute is
#' treated as one-dimensional.
#' @param intervals A list of minimum and maximum values. One for each
#' dimension of the hypercube.
#' @param threads An integer. In case `x` is a matrix of column vectors,
#' use this number of threads in parallel to evaluate.
#' @return A numeric. The interpolated value.
#' @export
#' @examples
#' # make a function which is known to be unsuitable for Chebyshev approximation
#' f <- function(x) sign(x)
#' # make a standard Chebyshev interpolation
#' ch <- ipol(f, dims = 50, method = "chebyshev")
#' # then do a truncated interpolation
#' val <- evalongrid(f, 50)
#' coef <- chebcoef(val)
#' # truncate the high frequencies
#' coef[-(1:10)] <- 0
#' # make a truncated approximation
#' tch <- Vectorize(function(x) chebeval(x, coef))
#' # make a lower degree also
#' ch2 <- ipol(f, dims = 10, method = "chebyshev")
#' # plot the functions
#' s <- seq(-1, 1, length.out = 400)
#' plot(s, ch(s), col = "red", type = "l")
#' lines(s, tch(s), col = "blue")
#' lines(s, f(s))
#' lines(s, ch2(s), col = "green")
chebeval <- function(x, coef, intervals = NULL, threads = getOption("chebpol.threads")) {
  if (is.integer(x) || is.logical(x)) {
    storage.mode(x) <- "double"
  }
  if (is.integer(coef) || is.logical(coef)) {
    storage.mode(coef) <- "double"
  }
  if (is.null(dim(coef))) {
    dim(coef) <- length(coef)
  }
  if (is.null(intervals)) {
    return(.Call(C_evalcheb, coef, x, threads, NULL))
  }
  if (is.numeric(intervals) && length(intervals) == 2L) {
    intervals <- list(intervals)
  }
  ispan <- vapply(intervals, \(i) 2 / diff(i), numeric(1L))
  mid <- vapply(intervals, mean, numeric(1L))
  .Call(C_evalcheb, coef, (x - mid) * ispan, threads, NULL)
}

# return a function which is a Chebyshev interpolation

chebappx.real <- function(val, intervals = NULL) {
  if (is.null(dim(val))) {
    # allow for one-dimensional
    dim(val) <- length(val)
  }
  K <- length(dim(val))
  # allow for vector, e.g. intervals = c(0, 1), put it inside list
  if (is.numeric(intervals) && length(intervals) == 2L) {
    intervals <- list(intervals)
  }
  cf <- chebcoef(val)

  if (is.null(intervals)) {
    # it's [-1, 1] intervals, so drop transformation
    fun <- vectorfun(
      function(x, threads = getOption("chebpol.threads")) .Call(C_evalcheb, cf, x, threads, NULL),
      arity = K,
      domain = rep(list(c(-1, 1)), K)
    )
  } else {
    # it's intervals, create mapping into [-1, 1]
    if (!is.list(intervals) || length(intervals) != K || any(lengths(intervals) != 2L)) {
      stop(sprintf("`intervals` must be a list of %d `c(lower, upper)` pairs, one per dimension of `val`.", K))
    }
    ispan <- vapply(intervals, \(x) 2 / diff(x), numeric(1L))
    mid <- vapply(intervals, mean, numeric(1L))
    imap <- function(x) (x - mid) * ispan

    fun <- vectorfun(
      function(x, threads = getOption("chebpol.threads")) .Call(C_evalcheb, cf, imap(x), threads, NULL),
      arity = K,
      domain = intervals
    )
  }
  rm(val)
  fun
}

chebappxf.real <- function(fun, dims, intervals = NULL, ...) {
  chebappx.real(evalongrid(fun, dims, intervals, ...), intervals)
}

# interpolate on a non-Chebyshev grid. This is useful if you for some reason
# do not have the function values on a Chebyshev-grid, but on some other grid.  It comes at a cost,
# The interpolation may not be very good compared to the Chebyshev-one.

# val are the function values on a grid, an array of appropriate dimension
# in the order of expand.grid()

# if grid is unspecified, it is assumed that it is on a Chebyshev grid in [-1, 1]
# If grid is specified, it is a list of vectors. The length of the list
# is the dimension of the grid. Each vector contains grid-points in increasing or decreasing order.
# val is assumed to be the function values on expand.grid(grid)
# The user-grid is mapped to Chebyshev knots in [-1, 1] by splinefun in pkg stats

chebappxg.real <- function(val, grid = NULL, mapdim = NULL) {
  # grid is a list of grid points. val is the values as in expand.grid(grid)
  # if grid is null it is assumed to be a chebyshev grid. The dimensions
  # must be present in val
  if (is.null(grid)) {
    return(chebappx.real(val))
  }
  if (is.null(dim(val))) {
    dim(val) <- length(val)
  }
  if (!is.list(grid) && length(grid) == length(val)) {
    grid <- list(grid)
  }
  if (prod(lengths(grid)) != length(val)) {
    stop(sprintf("`val` must have one value per grid point (%d), not %d.", prod(lengths(grid)), length(val)))
  }
  dim(val) <- lengths(grid)

  # ok, grid is something like list(c(...), c(...), c(...))
  # create piecewise linear functions which maps grid-points to chebyshev grids

  intervals <- lapply(grid, range)
  # create monotone splines with splinefun, method hyman
  gridmaps <- Map(stats::splinefun, grid, chebknots(dim(val)), MoreArgs = list(method = "hyman"))
  # map one dimension at a time, the spline functions are vectorized
  gridmap <- function(x) {
    if (!is.matrix(x)) {
      return(vapply(seq_along(gridmaps), \(i) gridmaps[[i]](x[i]), numeric(1L)))
    }
    for (i in seq_along(gridmaps)) {
      x[i, ] <- gridmaps[[i]](x[i, ])
    }
    x
  }
  ch <- chebappx.real(val)
  vectorfun(
    function(x, threads = getOption("chebpol.threads")) ch(gridmap(x), threads),
    arity = length(grid),
    domain = intervals
  )
}

chebappxgf.real <- function(fun, grid, ..., mapdim = NULL) {
  if (!is.list(grid)) {
    grid <- list(grid)
  }
  chebappxg.real(evalongrid(fun, ..., grid = grid), grid, mapdim)
}
