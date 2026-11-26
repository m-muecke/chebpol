# General grids, Floater-Hormann

ipol_fh <- function(val, grid = NULL, degree = 1L, ...) {
  if (!is.list(grid)) {
    grid <- list(grid)
  }
  grid <- lapply(grid, as.numeric)
  if (is.function(val)) {
    val <- evalongrid(val, grid = grid, ...)
  }
  grid_dims <- lengths(grid, use.names = FALSE)
  if (!is.array(val)) {
    dim(val) <- grid_dims
  } else if (!identical(unname(dim(val)), grid_dims)) {
    stop(sprintf(
      "`val` must have dimensions %s to match `grid`, not %s.",
      paste(grid_dims, collapse = " x "),
      paste(dim(val), collapse = " x ")
    ))
  }
  degree <- rep_len(as.integer(degree), length(grid))
  # calculate weights, formula 18 in Floater & Hormann, in C, parallelized
  weights <- .Call(c_fh_weights, grid, degree, getOption("chebpol.threads"))
  n_dim <- length(grid)
  fun <- function(x, threads = getOption("chebpol.threads")) {
    x <- as_points(x, n_dim)
    .Call(c_eval_fh, x, grid, val, weights, threads)
  }
  new_interpolant(crate(fun, val, grid, weights, n_dim), method = "fh", arity = n_dim, domain = lapply(grid, range))
}
