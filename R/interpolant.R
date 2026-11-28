#' Evaluate an interpolant in a point
#'
#' An interpolant is a function returned by [ipol()] which has prespecified values in some points,
#' and which fills in between with some reasonable values.
#'
#' Interpolants have class `chebpol_interpolant` and the attributes `method`, the interpolation
#' method, `arity`, the number of variables, and `domain`, a data frame with one row per variable
#' and the columns `variable`, `lower` and `upper`. The variables take their names from `grid`,
#' `dims`, `intervals` or the row names of `knots`, and are named `X1`, `X2`, ... otherwise.
#' Printing an interpolant shows these together with the defaults of its extra arguments, such as
#' `blend`.
#'
#' @name interpolant
#' @param x The argument of the function. A function of more than one variable takes a vector. `x`
#' can also be a matrix of column vectors.
#' @param threads The number of threads to use for evaluation. All interpolants created by chebpol
#' are parallelized. If given a matrix argument `x`, the vectors can be evaluated in parallel.
#' @param ... Other parameters. Currently used for simplex linear interpolants with the logical
#' argument `epol` which makes the interpolant extrapolate to points outside the domain. The
#' `"stalker"`, `"hstalker"`, `"oldstalker"`, `"multilinear"` and `"simplexlinear"` interpolants
#' have the argument `blend`, one of `"linear"`, `"cubic"`, `"sigmoid"`, `"parodic"` and `"square"`,
#' which chooses the blending function described in `vignette("stalker", package = "chebpol")`.
#' `"stalker"`, `"hstalker"` and `"simplexlinear"` also take `"mean"`. The default is `"cubic"` for
#' `"stalker"` and `"hstalker"`, and `"linear"` for the others.
#' @return A numeric. If more than one point was evaluated, a vector.
#' @examples
#' grid <- list(x = seq(0, 1, length.out = 10), y = seq(0, 1, length.out = 10))
#' val <- runif(100)
#' dim(val) <- c(10, 10)
#' ip <- ipol(val, grid = grid, method = "fh")
#' ip
#' ip(c(0.3, 0.8))
#' ip(matrix(runif(12), 2), threads = 2)
NULL

new_interpolant <- function(fun, method, arity, domain) {
  structure(
    fun,
    method = method,
    arity = arity,
    domain = as_domain(domain),
    class = c("chebpol_interpolant", "function")
  )
}

# One row per variable, named X1, X2, ... where no name was given
as_domain <- function(domain) {
  variable <- names(domain) %||% character(length(domain))
  unnamed <- is.na(variable) | !nzchar(variable)
  variable[unnamed] <- paste0("X", which(unnamed))
  data.frame(
    variable = variable,
    lower = vapply(domain, min, numeric(1L), USE.NAMES = FALSE),
    upper = vapply(domain, max, numeric(1L), USE.NAMES = FALSE)
  )
}

#' @export
print.chebpol_interpolant <- function(x, ...) {
  arity <- attr(x, "arity")
  cat(sprintf("<chebpol_interpolant> %s, %d variable%s\n", attr(x, "method"), arity, if (arity == 1L) "" else "s"))
  domain <- attr(x, "domain")
  lower <- vapply(domain$lower, format, character(1L))
  upper <- vapply(domain$upper, format, character(1L))
  bounds <- sprintf("[%s, %s]", lower, upper)
  cat("Domain: ", paste(bounds, collapse = " x "), "\n", sep = "")
  defaults <- formals(x)
  defaults <- defaults[setdiff(names(defaults), c("x", "threads"))]
  if (length(defaults) > 0L) {
    # evaluated in the interpolant, since `degree = default_degree` refers to a value there
    shown <- vapply(
      names(defaults),
      function(arg) {
        value <- eval(defaults[[arg]], environment(x))
        paste(arg, "=", paste(deparse(value), collapse = " "))
      },
      character(1L)
    )
    cat("Defaults: ", paste(shown, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
