slappx.real <- function(val, knots, ..., blend = NULL) {
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("method = \"simplexlinear\" requires the geometry package.")
  }
  if (is.function(val)) {
    val <- apply(knots, 2L, val, ...)
  }
  if (length(val) != ncol(knots)) {
    stop(sprintf("`val` must have one value per knot (%d), not %d.", ncol(knots), length(val)))
  }
  dtri <- t(geometry::delaunayn(t(knots), options = "Qt Pp"))
  # For fast evaluation: For each simplex, we precompute the LU-factorization
  # needed for transforming to barycentric coordinates. These are used for finding the right simplex.
  adata <- .Call(C_analyzesimplex, dtri, knots, getOption("chebpol.threads"))
  arity <- nrow(knots)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    epol = FALSE,
    blend = c("cubic", "linear", "sigmoid", "parodic", "square", "mean")
  ) {
    x <- as_points(x, arity)
    blend <- blend_types[[match.arg(blend)]]
    .Call(C_evalsl, x, knots, dtri, adata, val, epol, threads, blend)
  }
  new_interpolant(set_default_blend(fun, blend), method = "simplexlinear", arity = arity)
}
