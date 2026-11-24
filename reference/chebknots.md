# Create a Chebyshev-grid

Create a Chebyshev grid on a hypercube.

## Usage

``` r
chebknots(dims, intervals = NULL)
```

## Arguments

- dims:

  The number of grid-points in each dimension. For Chebyshev-polynomial
  of degree `dims-1`.

- intervals:

  A list of vectors of length 2. The lower and upper bounds of the
  hypercube.

## Value

A array of dimension `dims`. The Chebyshev grid-points.

## Details

If `intervals` is not provided, it is assumed that the domain of the
function in each dimension is `[-1, 1]`. Thus, standard Chebyshev knots
are produced. If `dims` is of length 1, `intervals` may be a vector of
length 2 rather than a list with a vector of length 2.

## Examples

``` r
# Standard knots for degree 3
chebknots(4)
#> [[1]]
#> [1]  0.9238795  0.3826834 -0.3826834 -0.9238795
#> 
# Knots in the interval [2,3] for degree 3
chebknots(4, interval = c(2, 3))
#> [[1]]
#> [1] 2.961940 2.691342 2.308658 2.038060
#> 
# Multivariate knots
chebknots(c(x = 3, y = 4, z = 3))
#> $x
#> [1]  0.8660254  0.0000000 -0.8660254
#> 
#> $y
#> [1]  0.9238795  0.3826834 -0.3826834 -0.9238795
#> 
#> $z
#> [1]  0.8660254  0.0000000 -0.8660254
#> 
# Multivariate grid
if (FALSE) { # \dontrun{
  expand.grid(chebknots(c(x = 3, y = 4, z = 5), list(c(1, 3), c(4, 6), c(800, 900))))
} # }
```
