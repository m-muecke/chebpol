# Evaluate a function on a grid

Evaluate a function on a Chebyshev grid, or on a user-specified grid.

## Usage

``` r
evalongrid(fun, dims, intervals = NULL, ..., grid = NULL)

evalongridV(fun, dims, intervals = NULL, ..., grid = NULL)
```

## Arguments

- fun:

  A multivariate real-valued function to be evaluated. Must be defined
  on the hypercube described by `intervals`.

- dims:

  A vector of integers. The number of grid points in each dimension.

- intervals:

  A list. Each entry is a vector of length 2 with the lower and upper
  end of the interval in each dimension.

- ...:

  Further arguments to `fun`.

- grid:

  Rather than specifying `dims` and `intervals` to get a Chebyshev grid,
  you may specify your own `grid` as a list of vectors whose Cartesian
  product will be the grid, as in `expand.grid(grid)`.

## Value

An array with the value of `fun` on each grid point. The `dim` attribute
has been appropriately set for the grid. If `fun` returns a vector, this
will be the first dimension of the returned array.

## Details

The function `fun` should be a `function(x, ...)`, where `length(x)`
equals `length(dims)` (or `length(grid)`).

If `grid` is provided, `fun` is evaluated on each point in the Cartesian
product of the vectors in `grid`.

If `intervals` is not provided, it is assumed that the domain of the
function is the hypercube `[-1, 1] x [-1, 1] x ... x [-1, 1]`. Thus, the
function is evaluated on a standard Chebyshev grid.

If `intervals` is provided, it should be a `list` with elements of
length 2, providing minimum and maximum for each dimension.

The grid itself may be produced by
`expand.grid(chebknots(dims, intervals))`, or `expand.grid(grid)`.

This function does the same as `apply(expand.grid(grid), 1, fun)`, but
it's faster and more memory-efficient for large grids because it does
not actually expand the grid.

The function `evalongridV()` is for vectorized functions, i.e. those
that can take a matrix of column vectors as argument. It's equivalent to
`fun(t(expand.grid(grid)))`.

## See also

Other Chebyshev tools:
[`chebcoef()`](https://m-muecke.github.io/chebpol/reference/chebcoef.md),
[`chebeval()`](https://m-muecke.github.io/chebpol/reference/chebeval.md),
[`chebknots()`](https://m-muecke.github.io/chebpol/reference/chebknots.md)

## Examples

``` r
f <- function(x) {
  a <- sum(x^2)
  ifelse(a == 0, 0, exp(-1 / a))
}
# Standard Chebyshev grid
evalongrid(f, dims = c(3, 5))
#>           [,1]       [,2]      [,3]       [,4]      [,5]
#> [1,] 0.5463972 0.40138578 0.2635971 0.40138578 0.5463972
#> [2,] 0.3310212 0.05533071 0.0000000 0.05533071 0.3310212
#> [3,] 0.5463972 0.40138578 0.2635971 0.40138578 0.5463972
# Then Chebyshev on [0,1] x [2,3]
evalongrid(f, dims = c(3, 5), intervals = list(c(0, 1), c(2, 3)))
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.9022755 0.8911381 0.8689764 0.8400496 0.8177098
#> [2,] 0.8959733 0.8832618 0.8574039 0.8224805 0.7945613
#> [3,] 0.8932498 0.8798220 0.8522416 0.8144193 0.7837026
# And on my own grid
grid <- list(sort(rnorm(3)), sort(rnorm(5)))
evalongrid(f, grid = grid)
#>           [,1]       [,2]        [,3]         [,4]         [,5]
#> [1,] 0.8178520 0.40472342 0.316158184 2.720667e-01 2.791339e-01
#> [2,] 0.7963855 0.14859337 0.030828349 4.788359e-03 7.193132e-03
#> [3,] 0.7891832 0.06025023 0.000222056 5.445309e-24 1.847445e-13
g <- ipol(f, grid = grid, method = "fh")
evalongridV(g, grid = grid, threads = 2)
#>           [,1]       [,2]        [,3]         [,4]         [,5]
#> [1,] 0.8178520 0.40472342 0.316158184 2.720667e-01 2.791339e-01
#> [2,] 0.7963855 0.14859337 0.030828349 4.788359e-03 7.193132e-03
#> [3,] 0.7891832 0.06025023 0.000222056 5.445309e-24 1.847445e-13
# vector valued function
f <- function(x) c(prod(x), sum(x^2))
evalongrid(f, grid = grid)
#> , , 1
#> 
#>          [,1]      [,2]       [,3]
#> [1,] 1.790345 0.8705785 -0.2189599
#> [2,] 4.973297 4.3922862  4.2237433
#> 
#> , , 2
#> 
#>           [,1]      [,2]        [,3]
#> [1,] 0.5120616 0.2489967 -0.06262536
#> [2,] 1.1055204 0.5245099  0.35596703
#> 
#> , , 3
#> 
#>           [,1]      [,2]        [,3]
#> [1,] 0.2859919 0.1390673 -0.03497693
#> [2,] 0.8684230 0.2874124  0.11886958
#> 
#> , , 4
#> 
#>            [,1]       [,2]         [,3]
#> [1,] 0.07446177 0.03620801 -0.009106707
#> [2,] 0.76822150 0.18721096  0.018668110
#> 
#> , , 5
#> 
#>            [,1]        [,2]       [,3]
#> [1,] -0.1315001 -0.06394364 0.01608252
#> [2,]  0.7836600  0.20264949 0.03410664
#> 
```
