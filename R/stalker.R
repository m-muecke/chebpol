# Precompute some values used in the stalker evaluation
stalkercontext <- function(dmin, dplus, r) {
  powmin <- Map(`^`, dmin, r)
  powplus <- Map(`^`, dplus, r)
  det <- Map(\(dm, dp, pm, pp) 1 / (dm * pp + dp * pm), dmin, dplus, powmin, powplus)
  list(det = det, pmin = powmin, pplus = powplus)
}

oldstalkerappx <- function(val, grid, r = 2, ..., blend = NULL) {
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
  gl <- prod(lengths(grid))
  if (length(val) != gl) {
    stop(sprintf("`val` must have one value per grid point (%d), not %d.", gl, length(val)))
  }

  r <- rep_len(as.numeric(r), length(grid))
  uniform <- vapply(grid, \(g) all(abs(diff(g, differences = 2L)) < 1e-8), logical(1L))
  # precompute distance between grid points
  dmin <- lapply(grid, \(gr) c(NaN, diff(gr)))
  dplus <- lapply(dmin, \(dd) c(dd[-1L], NaN))
  val <- as.numeric(val)
  stalker <- c(list(val = val, grid = grid), stalkercontext(dmin, dplus, r), uniform = list(uniform))

  if (!havegsl() && any(is.na(r) & !uniform)) {
    stop("`k = NA` on a non-uniform grid requires chebpol built with GSL.")
  }

  arity <- length(grid)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    degree = r,
    blend = c("linear", "cubic", "sigmoid", "parodic", "square")
  ) {
    x <- as_points(x, arity)
    blend <- blend_types[[match.arg(blend)]]
    if (!identical(degree, r)) {
      degree <- rep_len(as.numeric(degree), length(grid))
      if (!havegsl() && any(is.na(degree) & !uniform)) {
        stop("`degree = NA` on a non-uniform grid requires chebpol built with GSL.")
      }
      stalker <- c(list(val = val, grid = grid), stalkercontext(dmin, dplus, degree), uniform = list(uniform))
    }
    .Call(C_evalstalker, x, stalker, degree, blend, as.integer(threads))
  }
  new_interpolant(set_default_blend(fun, blend), method = "oldstalker", arity = arity, domain = lapply(grid, range))
}

hstalkerappx <- function(val, grid, ..., blend = NULL) {
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
  stalker <- .Call(C_makehyp, val, grid)
  arity <- length(grid)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    blend = c("cubic", "linear", "sigmoid", "parodic", "square", "mean")
  ) {
    x <- as_points(x, arity)
    blend <- blend_types[[match.arg(blend)]]
    .Call(C_evalhyp, x, stalker, blend, as.integer(threads))
  }
  new_interpolant(set_default_blend(fun, blend), method = "hstalker", arity = arity, domain = lapply(grid, range))
}

stalkerappx <- function(val, grid, ..., blend = NULL) {
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
  stalker <- .Call(C_makestalk, val, grid)
  arity <- length(grid)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    blend = c("cubic", "linear", "sigmoid", "parodic", "square", "mean")
  ) {
    x <- as_points(x, arity)
    blend <- blend_types[[match.arg(blend)]]
    .Call(C_evalstalk, x, stalker, blend, as.integer(threads))
  }
  new_interpolant(set_default_blend(fun, blend), method = "stalker", arity = arity, domain = lapply(grid, range))
}

havegsl <- function() .Call(C_havegsl)
