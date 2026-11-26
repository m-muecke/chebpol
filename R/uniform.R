# we can actually find the grid-maps for uniform grids.
# the Chebyshev knots are cos(pi * (j + 0.5) / n) for j = 0, ..., n - 1. These should
# map into the n grid points. These have distance 2 / (n - 1), starting in -1, ending in 1
# so they are -1 + 2 * j / (n - 1). After some manipulation, the function is:

ugm <- function(x, n) sinpi(0.5 * x * (1 - n) / n)

ucappx.real <- function(val, intervals = NULL) {
  if (is.null(dim(val))) {
    dim(val) <- length(val)
  }
  if (is.numeric(intervals) && length(intervals) == 2L) {
    intervals <- list(intervals)
  }
  dims <- dim(val)
  arity <- length(dims)
  cf <- chebcoef(val)
  # the per-dimension vectors recycle down the columns of x
  if (is.null(intervals)) {
    map <- crate(function(x) ugm(x, dims), dims)
  } else {
    assert_intervals(intervals, arity)
    # precompute interval mid points and inverse lengths
    mid <- vapply(intervals, mean, numeric(1L))
    ispan <- vapply(intervals, \(i) 2 / diff(i), numeric(1L))
    map <- crate(function(x) ugm(ispan * (x - mid), dims), ispan, mid, dims)
  }
  fun <- function(x, threads = getOption("chebpol.threads")) {
    x <- as_points(x, arity)
    .Call(c_eval_cheb, cf, map(x), threads, NULL)
  }
  domain <- stats::setNames(intervals %||% rep(list(c(-1, 1)), arity), names(intervals) %||% names(dims))
  new_interpolant(crate(fun, cf, map, arity), method = "uniform", arity = arity, domain = domain)
}

ucappxf.real <- function(fun, dims, intervals = NULL, ...) {
  if (is.numeric(intervals) && length(intervals) == 2L) {
    intervals <- list(intervals)
  }
  grid <- if (is.null(intervals)) {
    lapply(dims, \(d) seq(-1, 1, length.out = d))
  } else {
    assert_intervals(intervals, length(dims))
    Map(\(d, i) seq(i[1L], i[2L], length.out = d), dims, intervals)
  }
  val <- evalongrid(fun, ..., grid = grid)
  if (length(dim(val)) != length(dims)) {
    stop(sprintf("The function `val` must return a single number, not %d values.", dim(val)[1L]))
  }
  names(dim(val)) <- names(dims)
  ucappx.real(val, intervals)
}
