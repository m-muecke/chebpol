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
as_points <- function(x, arity, call = sys.call(sys.parent())) {
  if (is.integer(x) || is.logical(x)) {
    storage.mode(x) <- "double"
  }
  if (is.matrix(x)) {
    if (nrow(x) != arity) {
      msg <- sprintf("`x` must have %d rows, one per dimension, not %d.", arity, nrow(x))
      stop(errorCondition(msg, call = call))
    }
    return(x)
  }
  if (length(x) == arity) {
    return(x)
  }
  if (arity == 1L) {
    return(matrix(x, 1L))
  }
  if (length(x) %% arity != 0L) {
    msg <- sprintf("`x` must have length %d or be a matrix with %d rows, not length %d.", arity, arity, length(x))
    stop(errorCondition(msg, call = call))
  }
  n <- length(x) / arity
  if (n > 1L) {
    msg <- sprintf("Treating `x` of length %d as a %d x %d matrix, one point per column.", length(x), arity, n)
    warning(warningCondition(msg, call = call))
  }
  dim(x) <- c(arity, n)
  x
}

assert_intervals <- function(intervals, d, what = "val", call = sys.call(sys.parent())) {
  if (!is.list(intervals) || length(intervals) != d || any(lengths(intervals) != 2L)) {
    msg <- sprintf("`intervals` must be a list of %d `c(lower, upper)` pairs, one per dimension of `%s`.", d, what)
    stop(errorCondition(msg, call = call))
  }
}

is_flag <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

grid_knots <- function(grid) {
  knots <- t(expand.grid(lapply(grid, as.numeric)))
  rownames(knots) <- names(grid)
  knots
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs # nolint
}

blend_types <- c(linear = 0L, sigmoid = 1L, parodic = 2L, cubic = 3L, square = 4L, mean = 5L)

set_default_blend <- function(fun, blend) {
  if (is.null(blend)) {
    return(fun)
  }
  choices <- eval(formals(fun)$blend)
  blend <- match.arg(blend, choices)
  formals(fun)$blend <- c(blend, setdiff(choices, blend))
  fun
}
