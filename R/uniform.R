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
  ch <- chebappx.real(val)
  if (is.null(intervals)) {
    gridmap <- function(x) mapply(\(xi, d) ugm(xi, d), x, dims)
  } else {
    # precompute interval mid points and inverse lengths
    md <- lapply(intervals, mean)
    ispan <- lapply(intervals, \(i) 2 / diff(i))
    gridmap <- function(x) mapply(\(xi, mid, is, d) ugm(is * (xi - mid), d), x, md, ispan, dims)
  }
  gm <- function(x) {
    if (is.matrix(x)) apply(x, 2L, gridmap) else gridmap(x)
  }
  vectorfun(
    function(x, threads = getOption("chebpol.threads")) ch(gm(x), threads),
    arity = length(dims),
    domain = if (is.null(intervals)) rep(list(c(-1, 1)), length(dim(val))) else intervals
  )
}

ucappxf.real <- function(fun, dims, intervals = NULL, ...) {
  if (is.null(intervals)) {
    return(ucappx.real(evalongrid(fun, ..., grid = lapply(dims, \(d) seq(-1, 1, length.out = d)))))
  }
  if (is.numeric(intervals) && length(intervals) == 2L) {
    intervals <- list(intervals)
  }
  ucappx.real(
    evalongrid(
      fun,
      ...,
      grid = Map(\(d, i) seq(min(i), max(i), length.out = d), dims, intervals)
    ),
    intervals
  )
}
