# Precompute some values used in the stalker evaluation
stalker_context <- function(dmin, dplus, degree) {
  powmin <- Map(`^`, dmin, degree)
  powplus <- Map(`^`, dplus, degree)
  det <- Map(\(dm, dp, pm, pp) 1 / (dm * pp + dp * pm), dmin, dplus, powmin, powplus)
  list(det = det, pmin = powmin, pplus = powplus)
}

ipol_oldstalker <- function(val, grid, degree = 2, ..., blend = NULL) {
  if (!is.list(grid)) {
    grid <- list(grid)
  }
  grid <- lapply(grid, as.numeric)
  if (any(vapply(grid, is.unsorted, logical(1L)))) {
    if (!is.function(val)) {
      stop("`grid` must be increasing for method = \"oldstalker\".")
    }
    grid <- lapply(grid, sort)
  }
  if (is.function(val)) {
    val <- evalongrid(val, grid = grid, ...)
  }
  n <- prod(lengths(grid))
  if (length(val) != n) {
    stop(sprintf("`val` must have one value per grid point (%d), not %d.", n, length(val)))
  }

  default_degree <- rep_len(as.numeric(degree), length(grid))
  # relative to the extent of the grid, so that it does not depend on its scale
  uniform <- vapply(grid, \(g) all(abs(diff(g, differences = 2L)) < 1e-8 * diff(range(g))), logical(1L))
  # precompute distance between grid points
  dmin <- lapply(grid, \(gr) c(NaN, diff(gr)))
  dplus <- lapply(dmin, \(dd) c(dd[-1L], NaN))
  val <- as.numeric(val)
  stalker <- c(list(val = val, grid = grid), stalker_context(dmin, dplus, default_degree), uniform = list(uniform))

  if (!has_gsl() && any(is.na(default_degree) & !uniform)) {
    stop("`k = NA` on a non-uniform grid requires chebpol built with GSL.")
  }

  arity <- length(grid)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    degree = default_degree,
    blend = c("linear", "cubic", "sigmoid", "parodic", "square")
  ) {
    x <- as_points(x, arity)
    blend <- blend_types[[match.arg(blend)]]
    if (!identical(degree, default_degree)) {
      degree <- rep_len(as.numeric(degree), length(grid))
      if (!has_gsl() && any(is.na(degree) & !uniform)) {
        stop("`degree = NA` on a non-uniform grid requires chebpol built with GSL.")
      }
      stalker <- c(list(val = val, grid = grid), stalker_context(dmin, dplus, degree), uniform = list(uniform))
    }
    .Call(c_eval_oldstalker, x, stalker, degree, blend, as.integer(threads))
  }
  fun <- crate(fun, default_degree, grid, val, uniform, dmin, dplus, stalker, arity)
  new_interpolant(set_default_blend(fun, blend), method = "oldstalker", arity = arity, domain = lapply(grid, range))
}

ipol_hstalker <- function(val, grid, ..., blend = NULL) {
  if (!is.list(grid)) {
    grid <- list(grid)
  }
  grid <- lapply(grid, as.numeric)
  if (any(vapply(grid, is.unsorted, logical(1L)))) {
    if (!is.function(val)) {
      stop("`grid` must be increasing for method = \"hstalker\".")
    }
    grid <- lapply(grid, sort)
  }
  if (is.function(val)) {
    val <- evalongrid(val, grid = grid, ...)
  }
  stalker <- .Call(c_make_hstalker, val, grid)
  arity <- length(grid)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    blend = c("cubic", "linear", "sigmoid", "parodic", "square", "mean")
  ) {
    x <- as_points(x, arity)
    blend <- blend_types[[match.arg(blend)]]
    .Call(c_eval_hstalker, x, stalker, blend, as.integer(threads))
  }
  fun <- crate(fun, stalker, arity)
  new_interpolant(set_default_blend(fun, blend), method = "hstalker", arity = arity, domain = lapply(grid, range))
}

ipol_stalker <- function(val, grid, ..., blend = NULL) {
  if (!is.list(grid)) {
    grid <- list(grid)
  }
  grid <- lapply(grid, as.numeric)
  if (any(vapply(grid, is.unsorted, logical(1L)))) {
    if (!is.function(val)) {
      stop("`grid` must be increasing for method = \"stalker\".")
    }
    grid <- lapply(grid, sort)
  }
  if (is.function(val)) {
    val <- evalongrid(val, grid = grid, ...)
  }
  stalker <- .Call(c_make_stalker, val, grid)
  arity <- length(grid)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    blend = c("cubic", "linear", "sigmoid", "parodic", "square", "mean")
  ) {
    x <- as_points(x, arity)
    blend <- blend_types[[match.arg(blend)]]
    .Call(c_eval_stalker, x, stalker, blend, as.integer(threads))
  }
  fun <- crate(fun, stalker, arity)
  new_interpolant(set_default_blend(fun, blend), method = "stalker", arity = arity, domain = lapply(grid, range))
}
