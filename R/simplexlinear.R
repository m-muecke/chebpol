ipol_simplexlinear <- function(val, knots, ..., blend = NULL) {
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("method = \"simplexlinear\" requires the geometry package.")
  }
  if (is.null(dim(knots))) {
    dim(knots) <- c(1L, length(knots))
  }
  if (is.function(val)) {
    val <- apply(knots, 2L, val, ...)
  }
  if (length(val) != ncol(knots)) {
    stop(sprintf("`val` must have one value per knot (%d), not %d.", ncol(knots), length(val)))
  }
  if (is.integer(val) || is.logical(val)) {
    storage.mode(val) <- "double"
  }
  dtri <- t(geometry::delaunayn(t(knots), options = "Qt Pp"))
  # For fast evaluation: For each simplex, we precompute the LU-factorization
  # needed for transforming to barycentric coordinates. These are used for finding the right simplex.
  simplices <- .Call(c_make_simplexlinear, dtri, knots, getOption("chebpol.threads"))
  arity <- nrow(knots)
  fun <- function(
    x,
    threads = getOption("chebpol.threads"),
    epol = FALSE,
    blend = c("cubic", "linear", "sigmoid", "parodic", "square", "mean")
  ) {
    x <- as_points(x, arity)
    if (!is_flag(epol)) {
      stop("`epol` must be `TRUE` or `FALSE`.")
    }
    blend <- blend_types[[match.arg(blend)]]
    .Call(c_eval_simplexlinear, x, knots, dtri, simplices, val, epol, blend, threads)
  }
  fun <- crate(fun, knots, dtri, simplices, val, arity)
  new_interpolant(set_default_blend(fun, blend), method = "simplexlinear", arity = arity, domain = knots_domain(knots))
}
