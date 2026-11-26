ipol_crbf <- function(val, knots, rbase = 2, layers = 5, lambda = 0, ...) {
  if (is.null(dim(knots))) {
    dim(knots) <- c(1L, length(knots))
  }
  if (is.function(val)) {
    val <- apply(knots, 2L, val, ...)
  }
  model <- .Call(c_make_crbf, rbind(knots, val), layers, rbase, lambda)
  arity <- nrow(knots)
  fun <- function(x, threads = getOption("chebpol.threads")) {
    x <- as_points(x, arity)
    .Call(c_eval_crbf, x, model, threads)
  }
  new_interpolant(crate(fun, model, arity), method = "crbf", arity = arity, domain = knots_domain(knots))
}
