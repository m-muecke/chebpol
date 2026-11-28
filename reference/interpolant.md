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
  points outside the domain. The `"stalker"`, `"hstalker"`,
  `"oldstalker"`, `"multilinear"` and `"simplexlinear"` interpolants
  have the argument `blend`, one of `"linear"`, `"cubic"`, `"sigmoid"`,
  `"parodic"` and `"square"`, which chooses the blending function
  described in
  [`vignette("stalker", package = "chebpol")`](https://m-muecke.github.io/chebpol/articles/stalker.md).
  `"stalker"`, `"hstalker"` and `"simplexlinear"` also take `"mean"`.
  The default is `"cubic"` for `"stalker"` and `"hstalker"`, and
  `"linear"` for the others.

## Value

A numeric. If more than one point was evaluated, a vector.

## Details

Interpolants have class `chebpol_interpolant` and the attributes
`method`, the interpolation method, `arity`, the number of variables,
and `domain`, a data frame with one column per variable and the rows
`lower` and `upper`. The columns of `domain` take their names from
`grid`, `dims`, `intervals` or the row names of `knots`, and are named
`X1`, `X2`, ... otherwise. Printing an interpolant shows these together
with the defaults of its extra arguments, such as `blend`.

## Examples

``` r
grid <- list(x = seq(0, 1, length.out = 10), y = seq(0, 1, length.out = 10))
val <- runif(100)
dim(val) <- c(10, 10)
ip <- ipol(val, grid = grid, method = "fh")
ip
#> <chebpol_interpolant> fh, 2 variables
#> Domain: [0, 1] x [0, 1]
ip(c(0.3, 0.8))
#> [1] 0.6969905
ip(matrix(runif(12), 2), threads = 2)
#> [1] 0.2505217 1.2011162 1.2852535 0.6764973 1.0305336 1.3035894
```
