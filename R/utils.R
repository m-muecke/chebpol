`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs # nolint
}

crate <- function(.fn, ..., .parent = topenv(parent.frame())) {
  nn <- vapply(as.list(substitute(list(...)))[-1L], as.character, character(1L))
  environment(.fn) <- list2env(stats::setNames(list(...), nn), parent = .parent)
  .fn
}

knots_domain <- function(knots) {
  domain <- lapply(seq_len(nrow(knots)), \(i) range(knots[i, ]))
  names(domain) <- rownames(knots)
  domain
}

# `sys.call(-1L)` would name a helper when `as_points()` is forced lazily inside it
as_points <- function(x, n_dim, call = sys.call(sys.parent())) {
  if (is.integer(x) || is.logical(x)) {
    storage.mode(x) <- "double"
  }
  if (is.matrix(x)) {
    if (nrow(x) != n_dim) {
      msg <- sprintf("`x` must have %d rows, one per dimension, not %d.", n_dim, nrow(x))
      stop(errorCondition(msg, call = call))
    }
    return(x)
  }
  if (length(x) == n_dim) {
    return(x)
  }
  if (n_dim == 1L) {
    return(matrix(x, 1L))
  }
  if (length(x) %% n_dim != 0L) {
    msg <- sprintf("`x` must have length %d or be a matrix with %d rows, not length %d.", n_dim, n_dim, length(x))
    stop(errorCondition(msg, call = call))
  }
  n <- length(x) / n_dim
  if (n > 1L) {
    msg <- sprintf("Treating `x` of length %d as a %d x %d matrix, one point per column.", length(x), n_dim, n)
    warning(warningCondition(msg, call = call))
  }
  dim(x) <- c(n_dim, n)
  x
}

grid_knots <- function(grid) {
  knots <- t(expand.grid(lapply(grid, as.numeric)))
  rownames(knots) <- names(grid)
  knots
}

blend_types <- c(linear = 0L, sigmoid = 1L, parodic = 2L, cubic = 3L, square = 4L, mean = 5L)
