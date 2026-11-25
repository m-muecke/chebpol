rbf.alglib <- function(val, knots, rbase = 2, layers = 5, lambda = 0, ...) {
  if (is.null(dim(knots))) {
    dim(knots) <- c(1L, length(knots))
  }
  if (is.function(val)) {
    val <- apply(knots, 2L, val, ...)
  }
  model <- .Call(C_makerbf, rbind(knots, val), layers, rbase, lambda)
  arity <- nrow(knots)
  fun <- function(x, threads = getOption("chebpol.threads")) {
    x <- as_points(x, arity)
    .Call(C_evalrbf, model, x, threads)
  }
  new_interpolant(crate(fun, model, arity), method = "crbf", arity = arity, domain = knots_domain(knots))
}

#' Check whether chebpol has the ALGLIB library
#'
#' If ALGLIB was available at compile time, it can be used for compact
#' support radial basis function interpolations on scattered data.
#' This function checks whether ALGLIB is available.
#'
#' @return Returns `TRUE` if ALGLIB is available. Otherwise `FALSE`.
#' @export
havealglib <- function() .Call(C_havealglib)
