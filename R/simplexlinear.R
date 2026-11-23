slappx.real <- function(val, knots, ...) {
  if (is.function(val)) {
    val <- apply(knots, 2, val)
  }
  dtri <- t(geometry::delaunayn(t(knots), options = "Qt Pp"))
  # For fast evaluation: For each simplex, we precompute the LU-factorization
  # needed for transforming to barycentric coordinates. These are used for finding the right simplex.
  adata <- .Call(C_analyzesimplex, dtri, knots, getOption("chebpol.threads"))
  vectorfun(
    function(
      x,
      threads = getOption("chebpol.threads"),
      epol = FALSE,
      blend = c("cubic", "linear", "sigmoid", "parodic", "square", "mean")
    ) {
      blend <- switch(match.arg(blend), linear = 0L, sigmoid = 1L, parodic = 2L, cubic = 3L, square = 4L, mean = 5L)
      .Call(C_evalsl, x, knots, dtri, adata, val, epol, threads, as.integer(blend))
    },
    arity = nrow(knots)
  )
}
