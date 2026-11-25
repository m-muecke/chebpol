# Call after any change to the formals of `fun`, `formals<-` drops the class and attributes
new_interpolant <- function(fun, method, arity, domain = NULL) {
  structure(
    fun,
    method = method,
    arity = arity,
    domain = as_domain(domain),
    chebpol.version = utils::packageVersion("chebpol"),
    class = c("chebpol_interpolant", "function")
  )
}

crate <- function(.fn, ..., .parent = topenv(parent.frame())) {
  nn <- vapply(as.list(substitute(list(...)))[-1L], as.character, character(1L))
  environment(.fn) <- list2env(stats::setNames(list(...), nn), parent = .parent)
  .fn
}

# The ranges of the variables as a data frame, with the columns named X1, X2, ... where no name was given
as_domain <- function(domain) {
  if (length(domain) == 0L) {
    return(data.frame())
  }
  varnames <- names(domain) %||% character(length(domain))
  unnamed <- is.na(varnames) | varnames == ""
  varnames[unnamed] <- paste0("X", which(unnamed))
  names(domain) <- varnames
  as.data.frame(domain)
}

knots_domain <- function(knots) {
  domain <- lapply(seq_len(nrow(knots)), \(i) range(knots[i, ]))
  names(domain) <- rownames(knots)
  domain
}

# `sys.call(-1L)` would name a helper when `as_points()` is forced lazily inside it
as_points <- function(x, arity, call = sys.call(sys.parent())) {
  if (is.integer(x) || is.logical(x)) {
    storage.mode(x) <- "double"
  }
  if (is.matrix(x)) {
    if (nrow(x) != arity) {
      msg <- sprintf("`x` must have %d rows, one per dimension, not %d.", arity, nrow(x))
      stop(errorCondition(msg, call = call))
    }
    return(x)
  }
  if (length(x) == arity) {
    return(x)
  }
  if (arity == 1L) {
    return(matrix(x, 1L))
  }
  if (length(x) %% arity != 0L) {
    msg <- sprintf("`x` must have length %d or be a matrix with %d rows, not length %d.", arity, arity, length(x))
    stop(errorCondition(msg, call = call))
  }
  numvec <- length(x) / arity
  if (numvec > 1L) {
    msg <- sprintf("Treating `x` of length %d as a %d x %d matrix, one point per column.", length(x), arity, numvec)
    warning(warningCondition(msg, call = call))
  }
  dim(x) <- c(arity, numvec)
  x
}

grid_knots <- function(grid) {
  knots <- t(expand.grid(lapply(grid, as.numeric)))
  rownames(knots) <- names(grid)
  knots
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

blend_types <- c(linear = 0L, sigmoid = 1L, parodic = 2L, cubic = 3L, square = 4L, mean = 5L)

set_default_blend <- function(fun, blend) {
  if (is.null(blend)) {
    return(fun)
  }
  choices <- eval(formals(fun)$blend)
  blend <- match.arg(blend, choices)
  formals(fun)$blend <- c(blend, setdiff(choices, blend))
  fun
}

separgs <- function(fun, ...) {
  fun <- match.fun(fun)
  f <- function() {
    mc <- match.call()
    mc[[1L]] <- quote(list)
    arg <- eval.parent(mc)
    nn <- lengths(arg)
    if (all(nn == 1L)) {
      return(fun(as.numeric(arg)))
    }
    fun(do.call(rbind, arg))
  }
  args <- as.list(sys.call())[-(1:2)]
  formals(f) <- args
  f
}

#' Evaluate a function on a grid
#'
#' Evaluate a function on a Chebyshev grid, or on a user-specified grid.
#'
#' The function `fun` should be a `function(x, ...)`, where
#' `length(x)` equals `length(dims)` (or `length(grid)`).
#'
#' If `grid` is provided, `fun` is evaluated on each point in the
#' Cartesian product of the vectors in `grid`.
#'
#' If `intervals` is not provided, it is assumed that the domain of the
#' function is the hypercube `[-1, 1] x [-1, 1] x ... x [-1, 1]`. Thus, the
#' function is evaluated on a standard Chebyshev grid.
#'
#' If `intervals` is provided, it should be a `list` with elements of
#' length 2, providing minimum and maximum for each dimension.
#'
#' The grid itself may be produced by
#' `expand.grid(chebknots(dims, intervals))`, or
#' `expand.grid(grid)`.
#'
#' This function does the same as `apply(expand.grid(grid), 1, fun)`, but
#' it's faster and more memory-efficient for large grids because it does not
#' actually expand the grid.
#'
#' The function `evalongridV()` is for vectorized functions, i.e. those that
#' can take a matrix of column vectors as argument. It's equivalent to
#' `fun(t(expand.grid(grid)))`.
#'
#' @param fun Multivariate real-valued function to be evaluated. Must be
#' defined on the hypercube described by `intervals`.
#' @param dims A vector of integers. The number of grid-points in each
#' dimension.
#' @param intervals A list. Each entry is a vector of length 2 with the lower
#' and upper end of the interval in each dimension.
#' @param ... Further arguments to `fun`.
#' @param grid Rather than specifying `dims` and `intervals` to get a Chebyshev
#' grid, you may specify your own `grid` as a list of vectors whose
#' Cartesian product will be the grid, as in `expand.grid(grid)`.
#' @return An array with the value of `fun` on each grid point. The
#' `dim` attribute has been appropriately set for the grid. If `fun`
#' returns a vector, this will be the first dimension of the returned array.
#' @export
#' @examples
#' f <- function(x) {
#'   a <- sum(x^2)
#'   ifelse(a == 0, 0, exp(-1 / a))
#' }
#' # Standard Chebyshev grid
#' evalongrid(f, dims = c(3, 5))
#' # Then Chebyshev on [0,1] x [2,3]
#' evalongrid(f, dims = c(3, 5), intervals = list(c(0, 1), c(2, 3)))
#' # And on my own grid
#' grid <- list(sort(rnorm(3)), sort(rnorm(5)))
#' evalongrid(f, grid = grid)
#' g <- ipol(f, grid = grid, method = "fh")
#' evalongridV(g, grid = grid, threads = 2)
#' # vector valued function
#' f <- function(x) c(prod(x), sum(x^2))
#' evalongrid(f, grid = grid)
evalongrid <- function(fun, dims, intervals = NULL, ..., grid = NULL) {
  # do a call to stuff which doesn't really expand the grid
  if (is.numeric(grid)) {
    grid <- list(grid)
  }
  grid <- grid %||% chebknots(dims, intervals)
  mf <- match.fun(fun)
  .Call(C_evalongrid, \(x) mf(x, ...), lapply(grid, as.numeric))
}

#' @rdname evalongrid
#' @export
evalongridV <- function(fun, dims, intervals = NULL, ..., grid = NULL) {
  if (is.numeric(grid)) {
    grid <- list(grid)
  }
  grid <- grid %||% chebknots(dims, intervals)
  fun <- match.fun(fun)
  structure(fun(t(expand.grid(grid)), ...), dim = if (length(grid) > 1L) lengths(grid) else NULL)
}

#' Check whether chebpol uses FFTW
#'
#' It is possible to compile chebpol without FFTW. If this is done, it will
#' not be feasible to use high-degree Chebyshev polynomials. I.e. a 100x100x100
#' approximation will be possible, but not a one-dimensional 1000000. This
#' function checks whether chebpol uses FFTW.
#'
#' @return Returns `TRUE` if chebpol uses FFTW. Otherwise `FALSE`.
#' @export
havefftw <- function() .Call(C_havefftw)
