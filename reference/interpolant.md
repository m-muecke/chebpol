# Evaluate an interpolant in a point

An interpolant is a function returned by
[`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) which
has prespecified values in some points, and which fills in between with
some reasonable values.

## Arguments

- x:

  The argument of the function. A function of more than one variable
  takes a vector. `x` can also be a matrix of column vectors.

- threads:

  The number of threads to use for evaluation. All interpolants created
  by chebpol are parallelized. If given a matrix argument `x`, the
  vectors can be evaluated in parallel.

- ...:

  Other parameters. Currently used for simplex linear interpolants with
  the logical argument `epol` which makes the interpolant extrapolate to
  points outside the domain. The stalker spline has the argument
  `blend = c("linear", "cubic", "sigmoid")` where a blending function
  can be chosen as described in a vignette. The `"multilinear"`
  interpolant also has such a blending function.

## Value

A numeric. If more than one point was evaluated, a vector.

## Examples

``` r
grid <- list(x = seq(0, 1, length.out = 10), y = seq(0, 1, length.out = 10))
val <- runif(100)
dim(val) <- c(10, 10)
ip <- ipol(val, grid = grid, method = "fh")
ip(c(0.3, 0.8))
#> [1] 0.6969905
ip(matrix(runif(12), 2), threads = 2)
#> [1] 0.2505217 1.2011162 1.2852535 0.6764973 1.0305336 1.3035894
```
