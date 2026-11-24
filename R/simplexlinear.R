slappx.real <- function(val, knots, ...) {
  if (is.function(val)) {
    val <- apply(knots, 2L, val)
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
      blend <- blend_types[[match.arg(blend)]]
      .Call(C_evalsl, x, knots, dtri, adata, val, epol, threads, blend)
    },
    arity = nrow(knots)
  )
}
