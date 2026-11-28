#' Create an interpolating function
#'
#' Create an interpolating function from given values. Several interpolation methods are supported.
#'
#' `ipol()` is a wrapper around various interpolation methods in package chebpol. Which arguments to
#' specify depends on the method. The interpolation methods are described in `vignette("chebpol")`.
#'
#' The method `"chebyshev"` needs only the number of Chebyshev knots in each dimension. This is
#' either the `"dim"` attribute of the array `val`, or the `dims` argument if `val` is a function.
#' Also the intervals can be specified if different from `[-1, 1]`.
#'
#' The method `"uniform"` is similar to the `"chebyshev"`, but uniformly spaced knots are created.
#' The argument `intervals` generally goes with `dims` when other than standard intervals `[-1, 1]`
#' are used.
#'
#' The methods `"multilinear"`, `"fh"` (Floater-Hormann), `"stalker"`, `"hstalker"`, and `"general"`
#' need the argument `grid`. These are the methods which can use arbitrary Cartesian grids. The
#' stalker spline is described in `vignette("stalker")`. The Floater-Hormann method (`"fh"`) also
#' needs the `k` argument, the degree of the blending polynomials. It defaults to 4. It must be
#' between 0 and one less than the number of grid points in each dimension. Give either one value
#' for all dimensions or one value per dimension.
#'
#' The method `"polyharmonic"` needs the arguments `knots` and `k`. In addition it can take the
#' logical argument `normalize` for normalizing the knots to the unit hypercube. The default is
#' `NA`, which uses normalization if any of the knots are outside the unit hypercube. `normalize`
#' can also be a vector `c(a, b)`, or a matrix with columns `a` and `b` and one row per dimension,
#' in which case each coordinate is transformed as `a * (x - b)`. Also, a logical `nowarn` is
#' accepted. It is used to suppress a warning in case the system can't be solved exactly and a least
#' squares fallback method is used. Fitting solves a dense linear system with one row and column per
#' knot, so it becomes slow beyond a few thousand knots.
#'
#' The method `"simplexlinear"` needs the argument `knots`. It creates a Delaunay triangulation from
#' the knots, and does linear interpolation in each simplex by weighting the vertex values with the
#' barycentric coordinates.
#'
#' If knots are required, but the `grid` argument is given, knots are constructed as
#' `t(expand.grid(grid))`.
#'
#' There are also some usage examples and more in `vignette("chebpol")` and `vignette("chebusage")`.
#'
#' @param val An array or a function. Function values on a grid, or the function itself. If it is
#'   the values, the `"dim"` attribute must be appropriately set. If it is a function, it will be
#'   evaluated in the grid points.
#' @param dims An integer vector. The number of grid points in each dimension. Not needed if `val`
#'   is an array or `grid` is used.
#' @param intervals A list of numeric vectors of length 2. The lower and upper bound in each
#'   dimension. Not used if `grid` is specified.
#' @param grid A list. Each element is a vector of ordered grid points for a dimension. These need
#'   not be Chebyshev knots, nor evenly spaced.
#' @param knots A matrix. Each column is a point in an M-dimensional space.
#' @param k A numeric vector. Additional value, used with some methods.
#' @param method A string. The interpolation method to use.
#' @param ... Further arguments to the function, if `is.function(val)`. And some extra arguments for
#'   interpolant creation described in section Details.
#' @return A `function(x, threads = getOption("chebpol.threads"))` defined on a hypercube, an
#'   [interpolant] for the given function. The argument `x` can be a matrix of column vectors which
#'   are evaluated in parallel in a number of threads. The function yields values for arguments
#'   outside the hypercube as well, though it will typically be a poor approximation. `threads` is
#'   an integer specifying the number of parallel threads which should be used when evaluating a
#'   matrix of column vectors.
#' @export
#' @examples
#' # evenly spaced grid points
#' su <- seq(0, 1, length.out = 10)
#' # irregularly spaced grid points
#' s <- su^3
#' # create approximation on the irregularly spaced grid
#' ml1 <- ipol(exp, grid = list(s), method = "multilinear")
#' fh1 <- ipol(exp, grid = list(s), method = "fh")
#' # test it, since exp is convex, the linear approximation lies above
#' # the exp between the grid points
#' ml1(su) - exp(su)
#' fh1(su) - exp(su)
#'
#' # multi dimensional approximation
#' f <- function(x) 10 / (1 + 25 * mean(x^2))
#' # a 3-dimensional 10x10x10 grid, first and third coordinate are non-uniform
#' grid <- list(s, su, sort(1 - s))
#'
#' # make multilinear, Floater-Hormann, Chebyshev and polyharmonic spline.
#' ml2 <- ipol(f, grid = grid, method = "multilinear")
#' fh2 <- ipol(f, grid = grid, method = "fh")
#' hst <- ipol(f, grid = grid, method = "hstalker")
#' ch2 <- ipol(f, dims = c(10, 10, 10), intervals = list(0:1, 0:1, 0:1), method = "chebyshev")
#' knots <- matrix(runif(3 * 1000), 3)
#' ph2 <- ipol(f, knots = knots, k = 2, method = "polyharmonic")
#' # make 7 points in R3 to test them on
#' m <- matrix(runif(3 * 7), 3)
#' rbind(
#'   true = apply(m, 2, f),
#'   ml = ml2(m),
#'   fh = fh2(m),
#'   cheb = ch2(m),
#'   poly = ph2(m),
#'   hst = hst(m)
#' )
#' @examplesIf requireNamespace("geometry", quietly = TRUE)
#' # simplex linear interpolation needs the geometry package
#' sl2 <- ipol(f, knots = knots, method = "simplexlinear")
#' rbind(true = apply(m, 2, f), sl = sl2(m))
ipol <- function(
  val,
  dims = NULL,
  intervals = NULL,
  grid = NULL,
  knots = NULL,
  k = NULL,
  method = c(
    "chebyshev",
    "multilinear",
    "fh",
    "uniform",
    "general",
    "polyharmonic",
    "simplexlinear",
    "hstalker",
    "stalker",
    "oldstalker"
  ),
  ...
) {
  method <- match.arg(method)
  if (!is.function(val) && !is.numeric(val) && !is.logical(val)) {
    stop(sprintf("`val` must be numeric values or a function, not a %s.", class(val)[1L]))
  }
  if (is.integer(val) || is.logical(val)) {
    storage.mode(val) <- "double"
  }
  if (is.integer(knots)) {
    storage.mode(knots) <- "double"
  }
  # a single dimension may be given as a bare vector
  if (!is.null(grid) && !is.list(grid)) {
    grid <- list(grid)
  }
  if (!is.null(grid) && length(grid) == 0L) {
    stop("`grid` must have at least one dimension.")
  }
  if (!is.null(grid)) {
    grid <- lapply(grid, as.numeric)
  }
  if (!all(is.finite(unlist(grid, use.names = FALSE)))) {
    stop("`grid` must not contain missing or infinite values.")
  }
  switch(
    method,
    chebyshev = {
      if (is.function(val) && is.null(dims)) {
        stop("`dims` must be given when `val` is a function.")
      }
      if (!is.function(val) && is.null(dim(val)) && !is.null(dims)) {
        dim(val) <- dims
      }
      ipol_chebyshev(val, dims, intervals, ...)
    },
    multilinear = {
      assert_grid(grid, "multilinear", min_points = 2L)
      ipol_multilinear(val, grid, ...)
    },
    simplexlinear = {
      if (is.null(knots) && is.null(grid)) {
        stop("`knots` or `grid` must be given for method = \"simplexlinear\".")
      }
      knots <- knots %||% grid_knots(grid)
      ipol_simplexlinear(val, knots, ...)
    },
    fh = {
      assert_grid(grid, "fh")
      k <- k %||% pmin(4L, lengths(grid) - 1L)
      if (!length(k) %in% c(1L, length(grid))) {
        stop(sprintf("`k` must have length 1 or %d, not %d.", length(grid), length(k)))
      }
      k <- rep_len(k, length(grid))
      if (!is_whole(k) || any(k > lengths(grid) - 1L)) {
        stop("`k` must be whole numbers between 0 and one less than the number of grid points in each dimension.")
      }
      ipol_fh(val, grid, degree = k, ...)
    },
    uniform = {
      if (is.function(val) && is.null(dims)) {
        stop("`dims` must be given when `val` is a function.")
      }
      if (is.null(intervals) && !is.null(grid)) {
        intervals <- lapply(grid, range)
        warning("`intervals` not given, using the range of `grid`.")
      }
      if (!is.function(val) && is.null(dim(val)) && !is.null(dims)) {
        dim(val) <- dims
      }
      ipol_uniform(val, dims, intervals, ...)
    },
    general = {
      assert_grid(grid, "general")
      ipol_general(val, grid, ...)
    },
    polyharmonic = {
      if (is.null(knots) && is.null(grid)) {
        stop("`knots` or `grid` must be given for method = \"polyharmonic\".")
      }
      knots <- knots %||% grid_knots(grid)
      k <- k %||% 3
      ipol_polyharmonic(val, knots, k = k, ...)
    },
    oldstalker = {
      assert_grid(grid, "oldstalker", min_points = 2L)
      k <- k %||% 2
      ipol_oldstalker(val, grid, degree = k, ...)
    },
    hstalker = {
      assert_grid(grid, "hstalker", min_points = 2L)
      ipol_hstalker(val, grid, ...)
    },
    stalker = {
      assert_grid(grid, "stalker", min_points = 2L)
      ipol_stalker(val, grid, ...)
    }
  )
}
