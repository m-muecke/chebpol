# General grids, Floater-Hormann

fhappx.real <- function(val, grid = NULL, d = 1L, ...) {
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
  dd <- rep_len(as.integer(d), length(grid))
  # calculate weights, formula 18 in Floater & Hormann, in C, parallelized
  weights <- .Call(C_FHweights, grid, dd, getOption("chebpol.threads"))
  rm(dd, d)
  arity <- length(grid)
  new_interpolant(
    function(x, threads = getOption("chebpol.threads")) {
      x <- as_points(x, arity)
      .Call(C_FH, x, val, grid, weights, threads, NULL)
    },
    arity = arity,
    domain = lapply(grid, range)
  )
}
