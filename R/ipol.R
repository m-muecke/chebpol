#' Create interpolating function.
#'
#' Create an interpolating function from given values. Several interpolation methods are
#' supported.
#'
#' `ipol` is a wrapper around various interpolation methods in package \pkg{chebpol}.
#' Which arguments to specify depends on the method. The interpolation methods are described
#' in `vignette("chebpol", package = "chebpol")`.
#'
#' The method `"chebyshev"` needs only the number of Chebyshev knots in
#' each dimension. This is either the `"dim"` attribute of the array
#' `val`, or the `dims` argument if `val` is a function.  Also
#' the intervals can be specified if different from `[-1, 1]`.
#'
#' The method `"uniform"` is similar to the `"chebyshev"`, but
#' uniformly spaced knots are created.  The argument
#' `intervals` generally goes with `dims` when something else than
#' standard intervals `[-1, 1]` are used.
#'
#' The methods `"multilinear"`, `"fh"` (Floater-Hormann), `"stalker"`,
#' `"hstalker"`, and
#' `"general"` needs the argument `grid`.  These are the methods
#' which can use arbitrary Cartesian grids.
#' The stalker spline
#' is described in `vignette("stalker", package = "chebpol")`. The Floater-Hormann
#' method (`"fh"`) also needs the `k` argument, the degree of the blending
#' polynomials. It defaults to 4.
#'
#' The method `"polyharmonic"` needs the arguments `knots` and
#' `k`. In addition it can take the logical argument `normalize`
#' for normalizing the knots to the unit hypercube. The default is `NA`, which
#' uses normalization if any of the knots are outside the unit hypercube. Also, a logical
#' `nowarn` is accepted, it is used to suppress a warning in case the system can't
#' be solved exactly and a least squares fallback method is used.
#' Fitting solves a dense linear system with one row and column per knot, so it
#' becomes slow beyond a few thousand knots.
#'
#' The method `"simplexlinear"` needs the argument `knots`. It creates a
#' Delaunay triangulation from the knots, and does linear interpolation in each simplex
#' by weighting the vertex values with the barycentric coordinates.
#'
#' If knots are required, but the grid argument is given, knots are constructed as
#' `t(expand.grid(grid))`
#'
#' The `"crbf"` is the multilayer compact radial basis function
#' interpolation in ALGLIB (<http://www.alglib.net/interpolation/fastrbf.php>).
#' It is only available if ALGLIB was available at
#' compile time. It takes the extra arguments `"rbase"`, `"layers"`, and
#' `"lambda"`. These are discussed in the ALGLIB documentation.
#'
#' There are also some usage examples and more in `vignette("chebpol")` and `vignette("chebusage")`.
#'
#' @param val array or function. Function values on a grid, or the function
#' itself. If it is the values, the `"dim"`-attribute must be
#' appropriately set. If it is a function, it will be evaluated in the grid points.
#' @param dims integer vector. The number of grid points in each dimension. Not
#' needed if `val` is an array or `grid` is used.
#' @param intervals list of length 2 numeric vectors. The lower and upper bound
#' in each dimension. Not used if `grid` is specified.
#' @param grid list. Each element is a vector of ordered grid-points for a
#' dimension.  These need not be Chebyshev-knots, nor evenly spaced.
#' @param knots matrix. Each column is a point in an M-dimensional space.
#' @param k numeric. Additional value, used with some methods.
#' @param method character. The interpolation method to use.
#' @param ... Further arguments to the function, if `is.function(val)`. And some
#' extra arguments for interpolant creation described in section Details.
#' @return A `function(x, threads = getOption("chebpol.threads"))` defined on a hypercube, an [interpolant]
#' for the given function. The argument `x` can be a matrix of column
#' vectors which are evaluated in parallel in a number of threads.  The
#' function yields values for arguments outside the hypercube as well, though
#' it will typically be a poor approximation.  `threads` is an integer
#' specifying the number of parallel threads which should be used when
#' evaluating a matrix of column vectors.
#' @examples
#' # evenly spaced grid-points
#' su <- seq(0, 1, length.out = 10)
#' # irregularly spaced grid-points
#' s <- su^3
#' # create approximation on the irregularly spaced grid
#' ml1 <- ipol(exp, grid = list(s), method = "multilin")
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
#' ml2 <- ipol(f, grid = grid, method = "multilin")
#' fh2 <- ipol(f, grid = grid, method = "fh")
#' hst <- ipol(f, grid = grid, method = "hstalker")
#' ch2 <- ipol(f, dims = c(10, 10, 10), intervals = list(0:1, 0:1, 0:1), method = "cheb")
#' knots <- matrix(runif(3 * 1000), 3)
#' ph2 <- ipol(f, knots = knots, k = 2, method = "poly")
#' sl2 <- ipol(f, knots = knots, method = "simplexlinear")
#' # my alglib is a bit slow, so stick to 100 knots
#' if (havealglib()) {
#'   crb <- ipol(f, knots = knots[, 1:100], method = "crbf", rbase = 2, layers = 5, lambda = 0)
#' }
#' # make 7 points in R3 to test them on
#' m <- matrix(runif(3 * 7), 3)
#' rbind(
#'   true = apply(m, 2, f),
#'   ml = ml2(m),
#'   fh = fh2(m),
#'   cheb = ch2(m),
#'   poly = ph2(m),
#'   sl = sl2(m),
#'   hst = hst(m),
#'   crbf = if (havealglib()) crb(m) else NULL
#' )
#' @export
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
    "crbf",
    "oldstalker"
  ),
  ...
) {
  method <- match.arg(method)
  args <- list(...)
  switch(
    method,
    chebyshev = {
      if (is.function(val) && is.null(dims)) {
        stop("dims must be specified")
      }
      if (is.function(val)) {
        return(chebappxf.real(val, dims, intervals, ...))
      }
      if (is.null(dim(val)) && !is.null(dims)) {
        dim(val) <- dims
      }
      chebappx.real(val, intervals)
    },
    multilinear = {
      if (is.null(grid)) {
        stop("grid must be specified for multi linear interpolation")
      }
      if (!is.list(grid)) {
        grid <- list(grid)
      }
      if (any(lengths(grid) < 2L)) {
        stop("each grid dimension must have at least two points for multi linear interpolation")
      }
      grid <- lapply(grid, as.numeric)
      if (unsortedgrid(grid)) {
        stop("grid must be distinct ordered values")
      }
      blend <- args[["blend"]]
      if (is.null(blend)) {
        return(mlappx.real(val, grid, ...))
      }
      newarg <- args[!names(args) %in% "blend"]
      blenddef(do.call(mlappx.real, c(list(val, grid), newarg)), blend)
    },
    simplexlinear = {
      if (is.null(knots)) {
        if (!is.null(grid)) {
          knots <- t(expand.grid(grid))
        } else {
          stop("knots must be specified for simplex linear interpolation")
        }
      }
      blend <- args[["blend"]]
      if (is.null(blend)) {
        return(slappx.real(val, knots, ...))
      }
      newarg <- args[!names(args) %in% "blend"]
      blenddef(do.call(slappx.real, c(list(val, knots), newarg)), blend)
    },
    fh = {
      if (is.null(grid)) {
        stop("grid must be specified for Floater-Hormann interpolation")
      }
      if (!is.list(grid)) {
        grid <- list(grid)
      }
      if (unsortedgrid(grid)) {
        stop("grid must be distinct ordered values")
      }
      grid <- lapply(grid, as.numeric)
      k <- k %||% pmin(4L, lengths(grid) - 1L)
      fhappx.real(val, grid, d = k, ...)
    },
    uniform = {
      if (is.function(val) && is.null(dims)) {
        stop("Must specify dims for uniform intervals")
      }
      if (is.null(intervals) && !is.null(grid)) {
        intervals <- lapply(grid, range)
        warning("intervals constructed from ranges of grid-argument")
      }
      if (is.function(val)) {
        return(ucappxf.real(val, dims, intervals, ...))
      }
      if (is.null(dim(val)) && !is.null(dims)) {
        dim(val) <- dims
      }
      ucappx.real(val, intervals)
    },
    general = {
      if (is.null(grid)) {
        stop("grid must be specified for general interpolation")
      }
      if (!is.list(grid)) {
        grid <- list(grid)
      }
      if (unsortedgrid(grid)) {
        stop("grid must be distinct ordered values")
      }
      grid <- lapply(grid, as.numeric)
      if (is.function(val)) {
        return(chebappxgf.real(val, grid, ...))
      }
      chebappxg.real(val, grid)
    },
    polyharmonic = {
      if (is.null(knots)) {
        if (!is.null(grid)) {
          knots <- t(expand.grid(grid))
        } else {
          stop("Must specify knots for polyharmonic splines.")
        }
      }
      k <- k %||% 3
      # Now, ... may contain normalize and nowarn in addition to
      # arguments to val if it is a function
      # peel off normalize and nowarn
      normalize <- args[["normalize"]] %||% NA
      nowarn <- args[["nowarn"]] %||% FALSE
      newarg <- args[!names(args) %in% c("normalize", "nowarn")]
      do.call(polyh.real, c(list(val, knots, k, normalize, nowarn), newarg))
    },
    oldstalker = {
      if (is.null(grid)) {
        stop("grid must be specified for stalker interpolation")
      }
      if (!is.list(grid)) {
        grid <- list(grid)
      }
      if (unsortedgrid(grid)) {
        stop("grid must be distinct ordered values")
      }
      grid <- lapply(grid, as.numeric)
      k <- k %||% 2
      blend <- args[["blend"]]
      if (is.null(blend)) {
        return(oldstalkerappx(val, grid, r = k, ...))
      }
      newarg <- args[!names(args) %in% "blend"]
      blenddef(do.call(oldstalkerappx, c(list(val, grid, r = k), newarg)), blend)
    },
    hstalker = {
      if (is.null(grid)) {
        stop("grid must be specified for stalker interpolation")
      }
      if (!is.list(grid)) {
        grid <- list(grid)
      }
      if (unsortedgrid(grid)) {
        stop("grid must be distinct ordered values")
      }
      grid <- lapply(grid, as.numeric)
      blend <- args[["blend"]]
      if (is.null(blend)) {
        return(hstalkerappx(val, grid, ...))
      }
      newarg <- args[!names(args) %in% "blend"]
      blenddef(do.call(hstalkerappx, c(list(val, grid), newarg)), blend)
    },
    stalker = {
      if (is.null(grid)) {
        stop("grid must be specified for stalker interpolation")
      }
      if (!is.list(grid)) {
        grid <- list(grid)
      }
      if (unsortedgrid(grid)) {
        stop("grid must be distinct ordered values")
      }
      grid <- lapply(grid, as.numeric)
      blend <- args[["blend"]]
      if (is.null(blend)) {
        return(stalkerappx(val, grid, ...))
      }
      newarg <- args[!names(args) %in% "blend"]
      blenddef(do.call(stalkerappx, c(list(val, grid), newarg)), blend)
    },
    crbf = {
      if (is.null(knots)) {
        stop("Must specify knots for radial basis functions.")
      }
      if (length(args) == 0L) {
        # old interface
        if (is.null(k)) {
          rbase <- 2
          layers <- 5L
          lambda <- 0
        } else {
          rbase <- k[1L]
          layers <- k[2L]
          lambda <- k[3L]
        }
        return(rbf.alglib(val, knots, rbase, layers, lambda))
      }
      rbase <- args[["rbase"]] %||% 2
      layers <- args[["layers"]] %||% 5L
      lambda <- args[["lambda"]] %||% 0
      newarg <- args[!names(args) %in% c("rbase", "layers", "lambda")]
      do.call(rbf.alglib, c(list(val, knots, rbase, layers, lambda), newarg))
    },
    stop("Unknown interpolation method: ", method)
  )
}

unsortedgrid <- function(g) {
  any(vapply(g, \(s) is.unsorted(s, strictly = TRUE) && is.unsorted(-s, strictly = TRUE), logical(1L)))
}
