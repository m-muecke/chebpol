mlappx.real <- function(val, grid, ..., blend = NULL) {
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
  gl <- prod(lengths(grid))
  if (length(val) != gl) {
    stop(sprintf("`val` must have one value per grid point (%d), not %d.", gl, length(val)))
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
    .Call(C_evalmlip, grid, val, x, threads, blend)
  }
  new_interpolant(set_default_blend(fun, blend), arity = arity, domain = lapply(grid, range))
}
