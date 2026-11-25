# General grids, Floater-Hormann

fhappx.real <- function(val, grid = NULL, d = 1L, ...) {
  if (is.null(grid)) {
    stop("Must specify grid")
  }
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
      "dimensions of values (%s) must match the grid lengths (%s)",
      paste(dim(val), collapse = " x "),
      paste(grid_dims, collapse = " x ")
    ))
  }
  dd <- as.integer(d)
  dd <- rep(dd, length(grid) %/% length(d))
  # calculate weights, formula 18 in Floater & Hormann, in C, parallelized
  weights <- .Call(C_FHweights, grid, dd, getOption("chebpol.threads"))
  rm(dd, d)
  vectorfun(
    function(x, threads = getOption("chebpol.threads")) .Call(C_FH, x, val, grid, weights, threads, NULL),
    arity = length(grid),
    domain = lapply(grid, range)
  )
}
