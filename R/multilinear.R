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
  n_dim <- length(grid)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    blend = c("linear", "cubic", "sigmoid", "parodic", "square")
  ) {
    x <- as_points(x, n_dim)
    blend <- blend_types[[match.arg(blend)]]
    .Call(c_eval_multilinear, x, grid, val, blend, threads)
  }
  fun <- crate(fun, grid, val, n_dim)
  new_interpolant(set_default_blend(fun, blend), method = "multilinear", arity = n_dim, domain = lapply(grid, range))
}
