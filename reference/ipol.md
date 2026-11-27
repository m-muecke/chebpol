# Create an interpolating function

Create an interpolating function from given values. Several
interpolation methods are supported.

## Usage

``` r
ipol(
  val,
  dims = NULL,
  intervals = NULL,
  grid = NULL,
  knots = NULL,
  k = NULL,
  method = c("chebyshev", "multilinear", "fh", "uniform", "general", "polyharmonic",
    "simplexlinear", "hstalker", "stalker", "oldstalker"),
  ...
)
```

## Arguments

- val:

  array or function. Function values on a grid, or the function itself.
  If it is the values, the `"dim"`-attribute must be appropriately set.
  If it is a function, it will be evaluated in the grid points.

- dims:

  integer vector. The number of grid points in each dimension. Not
  needed if `val` is an array or `grid` is used.

- intervals:

  list of length 2 numeric vectors. The lower and upper bound in each
  dimension. Not used if `grid` is specified.

- grid:

  list. Each element is a vector of ordered grid-points for a dimension.
  These need not be Chebyshev-knots, nor evenly spaced.

- knots:

  matrix. Each column is a point in an M-dimensional space.

- k:

  numeric. Additional value, used with some methods.

- method:

  character. The interpolation method to use.

- ...:

  Further arguments to the function, if `is.function(val)`. And some
  extra arguments for interpolant creation described in section Details.

## Value

A `function(x, threads = getOption("chebpol.threads"))` defined on a
hypercube, an
[interpolant](https://m-muecke.github.io/chebpol/reference/interpolant.md)
for the given function. The argument `x` can be a matrix of column
vectors which are evaluated in parallel in a number of threads. The
function yields values for arguments outside the hypercube as well,
though it will typically be a poor approximation. `threads` is an
integer specifying the number of parallel threads which should be used
when evaluating a matrix of column vectors.

## Details

`ipol()` is a wrapper around various interpolation methods in package
chebpol. Which arguments to specify depends on the method. The
interpolation methods are described in
[`vignette("chebpol", package = "chebpol")`](https://m-muecke.github.io/chebpol/articles/chebpol.md).

The method `"chebyshev"` needs only the number of Chebyshev knots in
each dimension. This is either the `"dim"` attribute of the array `val`,
or the `dims` argument if `val` is a function. Also the intervals can be
specified if different from `[-1, 1]`.

The method `"uniform"` is similar to the `"chebyshev"`, but uniformly
spaced knots are created. The argument `intervals` generally goes with
`dims` when other than standard intervals `[-1, 1]` are used.

The methods `"multilinear"`, `"fh"` (Floater-Hormann), `"stalker"`,
`"hstalker"`, and `"general"` need the argument `grid`. These are the
methods which can use arbitrary Cartesian grids. The stalker spline is
described in
[`vignette("stalker", package = "chebpol")`](https://m-muecke.github.io/chebpol/articles/stalker.md).
The Floater-Hormann method (`"fh"`) also needs the `k` argument, the
degree of the blending polynomials. It defaults to 4. It must be between
0 and one less than the number of grid points in each dimension. Give
either one value for all dimensions or one value per dimension.

The method `"polyharmonic"` needs the arguments `knots` and `k`. In
addition it can take the logical argument `normalize` for normalizing
the knots to the unit hypercube. The default is `NA`, which uses
normalization if any of the knots are outside the unit hypercube.
`normalize` can also be a vector `c(a, b)`, or a matrix with columns `a`
and `b` and one row per dimension, in which case each coordinate is
transformed as `a * (x - b)`. Also, a logical `nowarn` is accepted. It
is used to suppress a warning in case the system can't be solved exactly
and a least squares fallback method is used. Fitting solves a dense
linear system with one row and column per knot, so it becomes slow
beyond a few thousand knots.

The method `"simplexlinear"` needs the argument `knots`. It creates a
Delaunay triangulation from the knots, and does linear interpolation in
each simplex by weighting the vertex values with the barycentric
coordinates.

If knots are required, but the `grid` argument is given, knots are
constructed as `t(expand.grid(grid))`.

There are also some usage examples and more in
[`vignette("chebpol")`](https://m-muecke.github.io/chebpol/articles/chebpol.md)
and
[`vignette("chebusage")`](https://m-muecke.github.io/chebpol/articles/chebusage.md).

## Examples

``` r
# evenly spaced grid-points
su <- seq(0, 1, length.out = 10)
# irregularly spaced grid-points
s <- su^3
# create approximation on the irregularly spaced grid
ml1 <- ipol(exp, grid = list(s), method = "multilinear")
fh1 <- ipol(exp, grid = list(s), method = "fh")
# test it, since exp is convex, the linear approximation lies above
# the exp between the grid points
ml1(su) - exp(su)
#>  [1] 0.0000000000 0.0007963441 0.0023666771 0.0036681447 0.0028930495
#>  [6] 0.0111173403 0.0064666457 0.0191999204 0.0246303860 0.0000000000
fh1(su) - exp(su)
#>  [1]  0.000000e+00  1.442078e-09 -2.352982e-08  7.217898e-08  1.287035e-07
#>  [6] -7.880319e-07 -8.460980e-07  3.720989e-06  1.076317e-05  0.000000e+00

# multi dimensional approximation
f <- function(x) 10 / (1 + 25 * mean(x^2))
# a 3-dimensional 10x10x10 grid, first and third coordinate are non-uniform
grid <- list(s, su, sort(1 - s))

# make multilinear, Floater-Hormann, Chebyshev and polyharmonic spline.
ml2 <- ipol(f, grid = grid, method = "multilinear")
fh2 <- ipol(f, grid = grid, method = "fh")
hst <- ipol(f, grid = grid, method = "hstalker")
ch2 <- ipol(f, dims = c(10, 10, 10), intervals = list(0:1, 0:1, 0:1), method = "chebyshev")
knots <- matrix(runif(3 * 1000), 3)
ph2 <- ipol(f, knots = knots, k = 2, method = "polyharmonic")
# make 7 points in R3 to test them on
m <- matrix(runif(3 * 7), 3)
rbind(
  true = apply(m, 2, f),
  ml = ml2(m),
  fh = fh2(m),
  cheb = ch2(m),
  poly = ph2(m),
  hst = hst(m)
)
#>          [,1]      [,2]     [,3]      [,4]      [,5]      [,6]      [,7]
#> true 1.798658 0.5777646 2.575931 0.7834249 0.8153359 0.6138369 1.0000857
#> ml   1.797827 0.5791828 2.595704 0.7825097 0.8157434 0.6173884 1.0001476
#> fh   1.798963 0.5782988 2.576696 0.7834285 0.8153104 0.6146149 1.0629866
#> cheb 1.798575 0.5777645 2.575911 0.7834264 0.8153354 0.6138369 1.0000685
#> poly 1.796145 0.5727733 2.575861 0.7841188 0.8153223 0.6139370 0.9992751
#> hst  1.795038 0.5817112 2.574213 0.7835482 0.8163775 0.6176764 1.0039954
# simplex linear interpolation needs the geometry package
sl2 <- ipol(f, knots = knots, method = "simplexlinear")
rbind(true = apply(m, 2, f), sl = sl2(m))
#>          [,1]      [,2]     [,3]      [,4]      [,5]      [,6]      [,7]
#> true 1.798658 0.5777646 2.575931 0.7834249 0.8153359 0.6138369 1.0000857
#> sl   1.805470        NA 2.574691 0.7865898 0.8160583 0.6129701 0.9950568
```
