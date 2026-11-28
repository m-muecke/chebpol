# Methods for creating multivariate interpolations on hypercubes

The package contains methods for creating multivariate/multidimensional
interpolations for real-valued functions on hypercubes. The methods
include classical Chebyshev interpolation, multilinear and
Floater-Hormann for arbitrary Cartesian product grids, and simplex
linear and polyharmonic spline for scattered multidimensional data.

## Details

The primary method of the package is
[`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) which
dispatches to some other method. All the generated
[interpolant](https://m-muecke.github.io/chebpol/reference/interpolant.md)s
accept as an argument a matrix of column vectors. The generated
functions also accept an argument
`threads = getOption("chebpol.threads")` to utilize more than one CPU if
a matrix of column vectors is evaluated. Unless the option
`chebpol.threads` is already set, it is initialized from the environment
variable `CHEBPOL_THREADS` upon loading of the package. It defaults to
`1`.

The interpolants are ordinary R-objects and can be saved with
[`save()`](https://rdrr.io/r/base/save.html) and loaded later with
[`load()`](https://rdrr.io/r/base/load.html) or serialized/unserialized
with other tools, just like any R-object. However, they call internal
functions of the package, which may change between versions, so an
interpolant saved with one version of chebpol may not work with another.
Create it again with
[`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) after
upgrading chebpol.

## Chebyshev

If we are free to evaluate the function to interpolate in arbitrary
points, we can use a Chebyshev interpolation. The classical one is
available with `ipol(..., method = "chebyshev")`.

## Uniform grids

There are several options if your function must be evaluated in a
uniform grid. There is the Floater-Hormann rational interpolation
available with `ipol(..., method = "fh")`. There is a transformed
Chebyshev variant `ipol(..., method = "uniform")`.

## Arbitrary grids

For grids which are not uniform, but still Cartesian products of
one-dimensional grids, there is the Floater-Hormann interpolation
`ipol(..., method = "fh")`, and a transformed Chebyshev variant
`ipol(..., method = "general")`, as well as a multilinear
`ipol(..., method = "multilinear")`. These methods work on uniform grids
as well. There is also a `method = "stalker"` which is a shape
preserving spline, and a variant `method = "hstalker"`. Both are
described in a vignette. These splines in one dimension attempt to
honour monotonicity and local extrema, in particular the `"hstalker"` is
positivity-preserving in the sense that it will preserve maxima and
minima.

## Scattered data

For scattered data, not necessarily organised as a Cartesian product
grid, there is a simplex linear interpolation available with
`ipol(..., method = "simplexlinear")`, and a polyharmonic spline with
`ipol(..., method = "polyharmonic")`.

## Support functions

There are also functions for producing Chebyshev grids
([`chebknots()`](https://m-muecke.github.io/chebpol/reference/chebknots.md))
as well as a function for evaluating a function on a grid
([`evalongrid()`](https://m-muecke.github.io/chebpol/reference/evalongrid.md)),
and a function for finding the Chebyshev coefficients
([`chebcoef()`](https://m-muecke.github.io/chebpol/reference/chebcoef.md)).

## See also

[`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md),
[interpolant](https://m-muecke.github.io/chebpol/reference/interpolant.md)

## Author

**Maintainer**: Maximilian Mücke <muecke.maximilian@gmail.com>
([ORCID](https://orcid.org/0009-0000-9432-9795))

Authors:

- Maximilian Mücke <muecke.maximilian@gmail.com>
  ([ORCID](https://orcid.org/0009-0000-9432-9795))

- Simen Gaure <Simen.Gaure@frisch.uio.no>
  ([ORCID](https://orcid.org/0000-0001-7251-8747)) \[copyright holder\]

## Examples

``` r
# make some function values on a 50x50x50 grid
dims <- c(x = 50, y = 50, z = 50)
f <- function(x) sqrt(1 + x[1]) * exp(x[2]) * sin(5 * x[3])^2
value <- evalongrid(f, dims)
# fit a Chebyshev approximation to it. Note that the value-array contains the grid-structure.
ch <- ipol(value, method = "chebyshev")
# To see the grid, use the chebknots function and expand.grid. Compare the
# values with the approximation in the first few grid points.
knots <- head(expand.grid(chebknots(dims)))
cbind(knots, value = head(as.numeric(value)), appx = ch(t(as.matrix(knots))))
#>           x         y         z    value     appx
#> 1 0.9995066 0.9995066 0.9995066 3.537865 3.537865
#> 2 0.9955620 0.9995066 0.9995066 3.534374 3.534374
#> 3 0.9876883 0.9995066 0.9995066 3.527394 3.527394
#> 4 0.9759168 0.9995066 0.9995066 3.516934 3.516934
#> 5 0.9602937 0.9995066 0.9995066 3.503002 3.503002
#> 6 0.9408808 0.9995066 0.9995066 3.485614 3.485614
# Make a Floater-Hormann approximation on a uniform grid as well
fh <- ipol(f, grid = lapply(dims, function(n) seq(-1, 1, length.out = n)), method = "fh", k = 5)
# evaluate in some random points in R3
m <- matrix(runif(15, -1, 1), 3)
rbind(true = apply(m, 2, f), cheb = ch(m), fh = fh(m))
#>           [,1]       [,2]       [,3]      [,4]      [,5]
#> true 0.1598858 0.03994209 0.05776106 0.4942973 0.3395190
#> cheb 0.1598858 0.03994154 0.05776073 0.4942993 0.3395198
#> fh   0.1598911 0.03994457 0.05775997 0.4942964 0.3395141
```
