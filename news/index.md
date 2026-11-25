# Changelog

## chebpol (development version)

### Breaking changes

- chebpol now requires R \>= 4.1.0.
- `chebappx()`, `chebappxf()`, `chebappxg()`, `chebappxgf()`,
  `fhappx()`, `mlappx()`, `polyh()`, `slappx()`, `ucappx()` and
  `ucappxf()`, deprecated since version 1.7, have been removed. Use
  [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md)
  instead.

### Bug fixes

- chebpol compiles again with current versions of R. FFTW is now found
  when it is installed outside the default include path.
- Integer input is now accepted everywhere a number is expected: the
  points passed to an interpolant, the values, grids and knots passed to
  [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md), the
  grid in
  [`evalongrid()`](https://m-muecke.github.io/chebpol/reference/evalongrid.md),
  and the arguments of
  [`chebcoef()`](https://m-muecke.github.io/chebpol/reference/chebcoef.md)
  and
  [`chebeval()`](https://m-muecke.github.io/chebpol/reference/chebeval.md).
  Previously most of these failed with “REAL() can only be applied to a
  ‘numeric’”.
- [`chebeval()`](https://m-muecke.github.io/chebpol/reference/chebeval.md)
  now works with `intervals` for a matrix of points, and for a single
  interval given as a vector `c(a, b)`. Before, only a single point in
  two or more dimensions was mapped correctly.
- Loading chebpol no longer overwrites the `chebpol.threads` option when
  it is already set, for example in `.Rprofile`. The `CHEBPOL_THREADS`
  environment variable is only used when the option is unset.
- [`evalongrid()`](https://m-muecke.github.io/chebpol/reference/evalongrid.md)
  now stops with an error when a grid dimension has no points.
  Previously it crashed R.
- [`evalongrid()`](https://m-muecke.github.io/chebpol/reference/evalongrid.md)
  and [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md)
  now accept functions that return integer or logical values, and stop
  with an error when the length of the returned value changes between
  grid points.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) now
  stops with an error when `blend` is not supported by the chosen
  method. Previously `blend = "mean"` was accepted for
  `method = "multilinear"` and `"oldstalker"` and silently treated as
  linear blending.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md)
  interpolants created with `method = "multilinear"` now return `NA` for
  `NA` or `NaN` input. Previously the evaluation never returned.
  [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) also
  stops with an error for a multilinear grid dimension with a single
  point, which was read out of bounds.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md)
  interpolants created with `method = "simplexlinear"` no longer return
  `NaN` with `blend = "square"` or `blend = "sigmoid"`. Square blending
  now uses the value at the vertex with the largest barycentric
  coordinate, and coordinates are clamped to `[0, 1]` before blending.
  `NA` or `NaN` input now gives `NA` for every blend.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "fh"` now stops with an error when the degree `k` is
  negative or not smaller than the number of grid points in a dimension.
  Before, the interpolant returned `NaN` between the grid points.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "fh"` now forwards extra arguments in `...` to the function
  being interpolated.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "uniform"` now accepts a single interval given as a vector
  `c(a, b)` together with values. It already did when given a function.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "simplexlinear"` now forwards extra arguments in `...` to
  the function being interpolated.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "polyharmonic"` accepts `normalize` as a vector `c(a, b)` or
  a two-column matrix again. Both forms failed with “the condition has
  length \> 1” on R \>= 4.2, and the knots were never transformed, so
  the fit did not match the evaluation.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md)
  interpolants created with `method = "polyharmonic"` now report the
  range of the knots in their `domain` attribute. With normalization it
  was the range after scaling, `[0, 1]`.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "polyharmonic"` no longer drops extra arguments in `...` for
  the function being interpolated when neither `normalize` nor `nowarn`
  is given.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "polyharmonic"` now computes the size of the knots-by-knots
  distance matrix correctly for 46,341 or more knots. Previously the
  size computation overflowed, giving “negative length vectors are not
  allowed” or, from 65,536 knots, possibly writing outside the allocated
  memory. Fitting that many knots still needs enough memory for the
  matrix, about 17 GB at 46,341 knots.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "polyharmonic"` now fits a spline with a single knot.
  Previously it failed with “invalid ‘nrow’ value”.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "polyharmonic"` now stops with an error when the number of
  values differs from the number of knots. Previously it warned about
  falling back to a least squares fit and suggested `normalize = TRUE`.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md)
  interpolants created with `method = "stalker"`, `"hstalker"` or
  `"simplexlinear"` now stop with an error when given a matrix with the
  wrong number of rows. Previously they returned values silently, and
  could read past the end of the input.
- [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md)
  interpolants created with `method = "oldstalker"` now report both
  ranks when given input of the wrong dimension.

### Documentation

- The vignettes are now R Markdown documents rendered to HTML instead of
  Sweave documents rendered to PDF.

## chebpol 2.0

CRAN release: 2018-09-15

- The stalker spline now works in N dimensions.
- A new hyperbolic stalker has also been introduced.
