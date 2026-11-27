is_flag <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

is_whole <- function(x, lower = 0) {
  is.numeric(x) && all(is.finite(x) & x == trunc(x) & x >= lower)
}

assert_intervals <- function(intervals, n_dim, what = "val", call = sys.call(sys.parent())) {
  if (!is.list(intervals) || length(intervals) != n_dim || any(lengths(intervals) != 2L)) {
    msg <- sprintf("`intervals` must be a list of %d `c(lower, upper)` pairs, one per dimension of `%s`.", n_dim, what)
    stop(errorCondition(msg, call = call))
  }
}

assert_grid <- function(grid, method, min_points = 1L, call = sys.call(sys.parent())) {
  if (is.null(grid)) {
    msg <- sprintf("`grid` must be given for method = \"%s\".", method)
    stop(errorCondition(msg, call = call))
  }
  if (any(lengths(grid) < min_points)) {
    msg <- sprintf(
      "Each `grid` dimension must have at least %s points for method = \"%s\".",
      c("one", "two")[min_points],
      method
    )
    stop(errorCondition(msg, call = call))
  }
  unsorted <- vapply(grid, \(g) is.unsorted(g, strictly = TRUE) && is.unsorted(-g, strictly = TRUE), logical(1L))
  if (any(unsorted)) {
    stop(errorCondition("Each `grid` dimension must be strictly increasing or strictly decreasing.", call = call))
  }
}
