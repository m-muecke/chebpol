# chebpol (development version)

## Breaking changes

* chebpol now requires R >= 4.1.0.
* `chebappx()`, `chebappxf()`, `chebappxg()`, `chebappxgf()`, `fhappx()`, `mlappx()`, `polyh()`, `slappx()`, `ucappx()` and `ucappxf()`, deprecated since version 1.7, have been removed. Use `ipol()` instead.

## Bug fixes

* chebpol compiles again with current versions of R. FFTW is now found when it is installed outside the default include path.
* Integer input is now accepted everywhere a number is expected: the points passed to an interpolant, the values, grids and knots passed to `ipol()`, the grid in `evalongrid()`, and the arguments of `chebcoef()` and `chebeval()`. Previously most of these failed with "REAL() can only be applied to a 'numeric'".
* `evalongrid()` and `ipol()` now accept functions that return integer or logical values, and stop with an error when the length of the returned value changes between grid points.
* `ipol()` now stops with an error when `blend` is not supported by the chosen method. Previously `blend = "mean"` was accepted for `method = "multilinear"` and `"oldstalker"` and silently treated as linear blending.
* `ipol()` interpolants created with `method = "multilinear"` now return `NA` for `NA` or `NaN` input. Previously the evaluation never returned. `ipol()` also stops with an error for a multilinear grid dimension with a single point, which was read out of bounds.
* `ipol()` interpolants created with `method = "simplexlinear"` no longer return `NaN` with `blend = "square"` or `blend = "sigmoid"`. Square blending now uses the value at the vertex with the largest barycentric coordinate, and coordinates are clamped to `[0, 1]` before blending. `NA` or `NaN` input now gives `NA` for every blend.
* `ipol()` with `method = "fh"` now forwards extra arguments in `...` to the function being interpolated.
* `ipol()` with `method = "simplexlinear"` now forwards extra arguments in `...` to the function being interpolated.
* `ipol()` with `method = "polyharmonic"` no longer drops extra arguments in `...` for the function being interpolated when neither `normalize` nor `nowarn` is given.
* `ipol()` interpolants created with `method = "stalker"`, `"hstalker"` or `"simplexlinear"` now stop with an error when given a matrix with the wrong number of rows. Previously they returned values silently, and could read past the end of the input.
* `ipol()` interpolants created with `method = "oldstalker"` now report both ranks when given input of the wrong dimension.

## Documentation

* The vignettes are now R Markdown documents rendered to HTML instead of Sweave documents rendered to PDF.

# chebpol 2.0

* The stalker spline now works in N dimensions.
* A new hyperbolic stalker has also been introduced.
