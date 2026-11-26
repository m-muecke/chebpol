ipol_multilinear <- function(val, grid, ..., blend = NULL) {
  if (is.numeric(grid)) {
    grid <- list(grid)
  }
  if (any(vapply(grid, is.unsorted, logical(1L)))) {
    if (!is.function(val)) {
      stop("`grid` must be increasing for method = \"multilinear\".")
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
  val <- as.numeric(val)
  arity <- length(grid)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    blend = c("linear", "cubic", "sigmoid", "parodic", "square")
  ) {
    x <- as_points(x, arity)
    blend <- blend_types[[match.arg(blend)]]
    .Call(c_eval_multilinear, grid, val, x, threads, blend)
  }
  fun <- crate(fun, grid, val, arity)
  new_interpolant(set_default_blend(fun, blend), method = "multilinear", arity = arity, domain = lapply(grid, range))
}
