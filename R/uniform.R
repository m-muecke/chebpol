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
    gm <- crate(function(x) ugm(x, dims), dims)
  } else {
    # precompute interval mid points and inverse lengths
    md <- vapply(intervals, mean, numeric(1L))
    ispan <- vapply(intervals, \(i) 2 / diff(i), numeric(1L))
    gm <- crate(function(x) ugm(ispan * (x - md), dims), ispan, md, dims)
  }
  fun <- function(x, threads = getOption("chebpol.threads")) {
    x <- as_points(x, arity)
    .Call(C_evalcheb, cf, gm(x), threads, NULL)
  }
  domain <- stats::setNames(intervals %||% rep(list(c(-1, 1)), arity), names(intervals) %||% names(dims))
  new_interpolant(crate(fun, cf, gm, arity), method = "uniform", arity = arity, domain = domain)
}

ucappxf.real <- function(fun, dims, intervals = NULL, ...) {
  if (is.numeric(intervals) && length(intervals) == 2L) {
    intervals <- list(intervals)
  }
  grid <- if (is.null(intervals)) {
    lapply(dims, \(d) seq(-1, 1, length.out = d))
  } else {
    Map(\(d, i) seq(i[1L], i[2L], length.out = d), dims, intervals)
  }
  val <- evalongrid(fun, ..., grid = grid)
  if (length(dim(val)) == length(dims)) {
    names(dim(val)) <- names(dims)
  }
  ucappx.real(val, intervals)
}
