ipol_crbf <- function(val, knots, rbase = 2, layers = 5, lambda = 0, ...) {
  if (is.null(dim(knots))) {
    dim(knots) <- c(1L, length(knots))
  }
  if (is.function(val)) {
    val <- apply(knots, 2L, val, ...)
  }
  model <- .Call(c_make_crbf, rbind(knots, val), layers, rbase, lambda)
  n_dim <- nrow(knots)
  fun <- function(x, threads = getOption("chebpol.threads")) {
    x <- as_points(x, n_dim)
    .Call(c_eval_crbf, x, model, threads)
  }
  new_interpolant(crate(fun, model, n_dim), method = "crbf", arity = n_dim, domain = knots_domain(knots))
}
