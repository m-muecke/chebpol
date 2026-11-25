# chebpol (development version)

## Breaking changes

* chebpol now requires R >= 4.1.0.
* `chebappx()`, `chebappxf()`, `chebappxg()`, `chebappxgf()`, `fhappx()`, `mlappx()`, `polyh()`, `slappx()`, `ucappx()` and `ucappxf()`, deprecated since version 1.7, have been removed. Use `ipol()` instead.
* geometry is now a suggested package instead of a required one, so installing chebpol no longer pulls in geometry and its dependencies. `ipol()` with `method = "simplexlinear"` needs it to build the triangulation and stops with an error when it is missing. Saved simplex linear interpolants can be evaluated without it.

## Minor improvements

* `ipol()` with `method = "polyharmonic"` now reports the warning about falling back to a least squares fit for the `ipol()` call. Previously it showed the source of an internal function, e.g. "In (function (val, knots, k = 2, ...".
* `ipol()` interpolants created with `method = "simplexlinear"` now record the range of their knots in the `domain` attribute, as polyharmonic interpolants do.
* `ipol()` interpolants now have class `chebpol_interpolant` and a `method` attribute. Printing one shows the method, the number of variables, the domain and the defaults of extra arguments such as `blend`, instead of the source code of the function.
* `ipol()` interpolants created with `method = "general"` or `"uniform"` evaluate a matrix of points about three times faster. The points are now mapped onto the Chebyshev grid one dimension at a time instead of one point at a time.
* `ipol()` interpolants take about half the time per call, which matters when they are called one point at a time. Errors about their input now name the interpolant, e.g. "Error in `ip()`" instead of "Error".

## Bug fixes

* `ipol()` interpolants no longer keep the extra arguments passed to the interpolated function, nor intermediate results such as the function values next to the Chebyshev coefficients. Previously an interpolant built with a large extra argument, e.g. `ipol(f, grid = grid, method = "fh", data = big)`, kept `big` in memory and in saved files.
* `ipol()` interpolants now name the columns of their `domain` attribute after the variables, taken from `grid`, `dims`, `intervals` or the row names of `knots`, and `X1`, `X2`, ... otherwise. Previously unnamed variables got names made from their bounds, such as `c.0..1.`, and names given in `dims` were dropped.
* chebpol compiles again with current versions of R. FFTW is now found when it is installed outside the default include path.
* Error messages are now consistent across interpolation methods and report the value that was given, e.g. "`x` must have 2 rows, not 3." instead of four differently worded messages for the same mistake.
* Integer input is now accepted everywhere a number is expected: the points passed to an interpolant, the values, grids and knots passed to `ipol()`, the grid in `evalongrid()`, and the arguments of `chebcoef()` and `chebeval()`. Previously most of these failed with "REAL() can only be applied to a 'numeric'".
* `chebeval()` now works with `intervals` for a matrix of points, and for a single interval given as a vector `c(a, b)`. Before, only a single point in two or more dimensions was mapped correctly.
* `chebcoef()` now stops with an error when `dct` is `NA` or not a single logical. Previously `dct = NA` was treated as `TRUE`.
* `chebeval()` now treats coefficients without a `dim` attribute, such as `c(1, 0.5, 0.25)` or `as.numeric(chebcoef(val))`, as one-dimensional. Previously it failed with "INTEGER() can only be applied to a 'integer', not a 'NULL'".
* `chebeval()` now errors when `coef` is empty. Previously it could crash R.
* Loading chebpol no longer overwrites the `chebpol.threads` option when it is already set, for example in `.Rprofile`. The `CHEBPOL_THREADS` environment variable is only used when the option is unset.
* `evalongrid()` and `ipol()` no longer overwrite an argument vector that the function has kept. Previously one vector was reused for every grid point, so a function that kept its argument, for example to record the points it was called with, saw every copy change to the last grid point.
* `evalongrid()` now stops with an error when a grid dimension has no points. Previously it crashed R.
* `evalongrid()` and `ipol()` now stop with an error when `grid` has no dimensions, e.g. `grid = list()`. Previously `ipol()` returned an interpolant without variables for several methods, and other methods failed with "'dim' cannot be of length 0".
* `evalongrid()` and `ipol()` now accept functions that return integer or logical values, and stop with an error when the length of the returned value changes between grid points.
* `ipol()` interpolants created with an explicit `blend` now keep their `arity`, `domain` and `chebpol.version` attributes. Previously setting the default blend dropped them.
* `ipol()` now stops with an error when `val` is neither numeric values nor a function, such as a function name given as a string. Previously such input failed with unrelated errors, and a list of values was accepted by some methods.
* `ipol()` now stops with an error when `blend` is not supported by the chosen method. Previously `blend = "mean"` was accepted for `method = "multilinear"` and `"oldstalker"` and silently treated as linear blending.
* `ipol()` interpolants created with `method = "multilinear"` now return `NA` for `NA` or `NaN` input. Previously the evaluation never returned. `ipol()` also stops with an error for a multilinear grid dimension with a single point, which was read out of bounds.
* `ipol()` interpolants created with `method = "simplexlinear"` no longer return `NaN` with `blend = "square"` or `blend = "sigmoid"`. Square blending now uses the value at the vertex with the largest barycentric coordinate, and coordinates are clamped to `[0, 1]` before blending. `NA` or `NaN` input now gives `NA` for every blend.
* `ipol()` with `method = "fh"` now stops with an error when the degree `k` is negative or not smaller than the number of grid points in a dimension. Before, the interpolant returned `NaN` between the grid points.
* `ipol()` with `method = "fh"` now stops with an error when `k` has neither length 1 nor one value per grid dimension. Previously a shorter `k` was silently repeated, or failed with "Length of grid (3) should equal length of d (2)".
* `ipol()` with `method = "fh"` now forwards extra arguments in `...` to the function being interpolated.
* `ipol()` with `method = "fh"` now stops with an error when the dimensions of the values do not match the grid lengths. Previously the interpolant read past the end of the grid and returned meaningless values.
* `ipol()` with `method = "uniform"` now accepts a one-dimensional `grid` given as a bare vector, as the other grid methods do. Previously the interpolant returned `NaN`s.
* `ipol()` with `method = "uniform"` now handles a reversed interval such as `c(1, 0)` when given a function. Previously the interpolant was silently wrong.
* `ipol()` with `method = "uniform"` now accepts a single interval given as a vector `c(a, b)` together with values. It already did when given a function.
* `ipol()` with `method = "simplexlinear"` now forwards extra arguments in `...` to the function being interpolated.
* `ipol()` with `method = "simplexlinear"` now stops with an error when the number of values differs from the number of knots. Previously the interpolant read past the end of the values and returned meaningless values or `NA`.
* `ipol()` with `method = "polyharmonic"` accepts `normalize` as a vector `c(a, b)` or a two-column matrix again. Both forms failed with "the condition has length > 1" on R >= 4.2, and the knots were never transformed, so the fit did not match the evaluation.
* `ipol()` interpolants created with `method = "polyharmonic"` now report the range of the knots in their `domain` attribute. With normalization it was the range after scaling, `[0, 1]`.
* `ipol()` with `method = "polyharmonic"` no longer drops extra arguments in `...` for the function being interpolated when neither `normalize` nor `nowarn` is given.
* `ipol()` with `method = "polyharmonic"` now computes the size of the knots-by-knots distance matrix correctly for 46,341 or more knots. Previously the size computation overflowed, giving "negative length vectors are not allowed" or, from 65,536 knots, possibly writing outside the allocated memory. Fitting that many knots still needs enough memory for the matrix, about 17 GB at 46,341 knots.
* `ipol()` with `method = "polyharmonic"` now fits knots on a line or plane, such as a single knot or knots sharing a coordinate, and knots with a large offset. Previously these returned `NaN` or fell back to an inaccurate least squares fit with a warning.
* `ipol()` with `method = "polyharmonic"` now fits a spline with a single knot. Previously it failed with "invalid 'nrow' value".
* `ipol()` with `method = "polyharmonic"` now stops with an error when the number of values differs from the number of knots. Previously it warned about falling back to a least squares fit and suggested `normalize = TRUE`.
* `ipol()` with `method = "stalker"`, `"hstalker"` or `"oldstalker"` now stops with an error when a grid dimension has a single point. Previously the spline read past the end of the values and the grid.
* `ipol()` interpolants created with `method = "hstalker"` no longer return `NaN` at a grid point with `blend = "mean"` when two neighboring values are equal.
* `ipol()` with `method = "hstalker"` no longer depends on the scale of the values or of the grid. Previously small differences between neighboring values were treated as constant or linear, and could produce a narrow spike between grid points.
* `ipol()` interpolants created with `method = "stalker"`, `"hstalker"` or `"simplexlinear"` now stop with an error when given a matrix with the wrong number of rows. Previously they returned values silently, and could read past the end of the input.
* `ipol()` interpolants created with `method = "oldstalker"` now report both ranks when given input of the wrong dimension.

## Documentation

* The vignettes are now R Markdown documents rendered to HTML instead of Sweave documents rendered to PDF.

# chebpol 2.0

* The stalker spline now works in N dimensions.
* A new hyperbolic stalker has also been introduced.
