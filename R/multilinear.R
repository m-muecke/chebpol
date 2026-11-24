mlappx.real <- function(val, grid, ...) {
  if (is.numeric(grid)) {
    grid <- list(grid)
  }
  if (any(vapply(grid, is.unsorted, logical(1L)))) {
    if (!is.function(val)) {
      stop("Grid points must be ordered in increasing order")
    }
    grid <- lapply(grid, sort)
  }
  if (is.function(val)) {
    val <- evalongrid(val, grid = grid, ...)
  }
  gl <- prod(lengths(grid))
  if (length(val) != gl) {
    stop("length of values ", length(val), " do not match size of grid ", gl)
  }
  val <- as.numeric(val)
  vectorfun(
    function(x, threads = getOption("chebpol.threads"), blend = c("linear", "cubic", "sigmoid", "parodic", "square")) {
      blend <- blend_types[[match.arg(blend)]]
      .Call(C_evalmlip, grid, val, x, threads, blend)
    },
    arity = length(grid),
    domain = lapply(grid, range)
  )
}
