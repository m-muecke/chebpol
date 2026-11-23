#' Evaluate an interpolant in a point
#'
#' An interpolant is a function returned by [ipol()] which has prespecified values in some
#' points, and which fills in between with some reasonable values.
#'
#' @name interpolant
#' @param x The argument of the function. A function of more then one variable takes a
#' vector. `x` can also be a matrix of column vectors.
#' @param threads The number of threads to use for evaluation. All  interpolants created by
#' \pkg{chebpol} are parallelized. If given a matrix argument `x`, the vectors can
#' be evaluated in parallel.
#' @param ... Other parameters. Currently used for simplex linear interpolants with the logical argument
#' `epol` which makes the interpolant extrapolate to points outside the domain.
#' The stalker spline has the argument `blend=c("linear","cubic","sigmoid")` where a
#' blending function can be chosen as described in a vignette. The `"multilinear"` interpolant also
#' has such a blending function.
#' @return A numeric. If more than one point was evaluated, a vector.
#' @examples
#' grid <- list(x = seq(0, 1, length.out = 10), y = seq(0, 1, length.out = 10))
#' val <- runif(100)
#' dim(val) <- c(10, 10)
#' ip <- ipol(val, grid = grid, method = "fh")
#' ip(c(0.3, 0.8))
#' ip(matrix(runif(12), 2), threads = 2)
NULL
