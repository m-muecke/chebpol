# chebpol (development version)

## Breaking changes

* chebpol now requires R >= 4.1.0.
* `chebappx()`, `chebappxf()`, `chebappxg()`, `chebappxgf()`, `fhappx()`, `mlappx()`, `polyh()`, `slappx()`, `ucappx()` and `ucappxf()`, deprecated since version 1.7, have been removed. Use `ipol()` instead.
* `ipol()` interpolants no longer have a `chebpol.version` attribute. Nothing in chebpol used it, and looking up the version took most of the time needed to create a small interpolant, which is now about four times faster.
* `ipol()` interpolants saved with an earlier version of chebpol can no longer be evaluated, because the internal functions they call have been renamed. Create them again with `ipol()`.
* The C routines `evalcheb`, `evalmlip` and `evalongrid`, registered for other packages to call through `R_GetCCallable()`, have been removed. They were undocumented, had no header file, and no known package used them.
* geometry is now a suggested package instead of a required one, so installing chebpol no longer pulls in geometry and its dependencies. `ipol()` with `method = "simplexlinear"` needs it to build the triangulation and stops with an error when it is missing. Saved simplex linear interpolants can be evaluated without it.
* `ipol()` interpolants created with `method = "simplexlinear"` now use `blend = "linear"` by default instead of `"cubic"`, so they interpolate linearly in each simplex, as documented.
* `ipol()` no longer has `method = "crbf"`. It needed ALGLIB, which chebpol has not supported since version 2.1, so it always stopped with an error.

## Lifecycle changes

* `havefftw()` is deprecated in favor of `has_fftw()`. `havealglib()` is deprecated and always returns `FALSE`. The new `has_gsl()` reports whether chebpol was built with GSL, which the stalker methods need to solve for the degree on a non-uniform grid.

## Minor improvements

* Attaching chebpol no longer prints a message when it was built without FFTW. Use `has_fftw()` to check.
* chebpol can be built without FFTW or GSL by setting the environment variable `CHEBPOL_DISABLE_FFTW` or `CHEBPOL_DISABLE_GSL` before installing, as an alternative to `--configure-args`.
* `ipol()` with `method = "polyharmonic"` now reports the warning about falling back to a least squares fit for the `ipol()` call. Previously it showed the source of an internal function, e.g. "In (function (val, knots, k = 2, ...".
* `ipol()` interpolants created with `method = "multilinear"`, `"stalker"` or `"hstalker"` evaluate faster in several dimensions with a blend other than `"linear"`, e.g. about four times faster for `"multilinear"` and twice as fast for the stalker methods with `blend = "sigmoid"` in five dimensions. The blending weights are now computed once per dimension instead of once per corner of the grid cell.
* `ipol()` interpolants created with `method = "simplexlinear"` evaluate much faster with many knots, e.g. about 200 times faster with 10,000 knots in two or three dimensions. The simplex holding a point is now looked up in a grid instead of by scanning all simplices, so the time per point no longer grows with the number of knots in up to three dimensions. The gain is smaller in higher dimensions.
* `ipol()` interpolants created with `method = "simplexlinear"` now record the range of their knots in the `domain` attribute, as polyharmonic interpolants do.
* `ipol()` interpolants now have class `chebpol_interpolant` and a `method` attribute. Printing one shows the method, the number of variables, the domain and the defaults of extra arguments such as `blend`, instead of the source code of the function.
* `ipol()` interpolants created with `method = "general"` or `"uniform"` evaluate a matrix of points about three times faster. The points are now mapped onto the Chebyshev grid one dimension at a time instead of one point at a time.
* `ipol()` interpolants take about half the time per call, which matters when they are called one point at a time. Errors about their input now name the interpolant, e.g. "Error in `ip()`" instead of "Error".

## Bug fixes

* `ipol()` interpolants no longer keep the extra arguments passed to the interpolated function, nor intermediate results such as the function values next to the Chebyshev coefficients. Previously an interpolant built with a large extra argument, e.g. `ipol(f, grid = grid, method = "fh", data = big)`, kept `big` in memory and in saved files.
* `ipol()` interpolants now have a `domain` attribute with one row per variable and the columns `name`, `lower` and `upper`. The variables are named after `grid`, `dims`, `intervals` or the row names of `knots`, and `X1`, `X2`, ... otherwise. Previously it had one column per variable, with names made from the bounds for unnamed variables, such as `c.0..1.`, and names given in `dims` were dropped.
* chebpol compiles again with current versions of R. FFTW is now found when it is installed outside the default include path.
* Error messages are now consistent across interpolation methods and report the value that was given, e.g. "`x` must have 2 rows, not 3." instead of four differently worded messages for the same mistake.
* Integer input is now accepted everywhere a number is expected: the points passed to an interpolant, the values, grids and knots passed to `ipol()`, the grid in `evalongrid()`, and the arguments of `chebcoef()` and `chebeval()`. Previously most of these failed with "REAL() can only be applied to a 'numeric'".
* `ipol()` with `method = "fh"`, `"stalker"` or `"hstalker"` now accepts logical values, as the other methods do. Previously these failed with "REAL() can only be applied to a 'numeric', not a 'logical'".
* `chebeval()` now works with `intervals` for a matrix of points, and for a single interval given as a vector `c(a, b)`. Before, only a single point in two or more dimensions was mapped correctly.
* `chebcoef()` now stops with an error when `dct` is `NA` or not a single logical. Previously `dct = NA` was treated as `TRUE`.
* `chebeval()` and `ipol()` with `method = "uniform"` now stop with an error when `intervals` is not a list with one `c(lower, upper)` pair per dimension, as `method = "chebyshev"` does. Previously a list of the wrong length was silently recycled, which mapped the points wrongly or gave an interpolant that could not be evaluated.
* `chebeval()` now treats coefficients without a `dim` attribute, such as `c(1, 0.5, 0.25)` or `as.numeric(chebcoef(val))`, as one-dimensional. Previously it failed with "INTEGER() can only be applied to a 'integer', not a 'NULL'".
* `chebeval()` now errors when `coef` is empty. Previously it could crash R.
* `chebknots()`, `evalongrid()` and `ipol()` with `method = "chebyshev"` now stop with an error when `dims` is not whole numbers, or when `intervals` is not a list with one `c(lower, upper)` pair per dimension. Previously `dims = 2.5` gave two knots that are not Chebyshev knots, so `ipol()` silently returned a wrong interpolant, and a list of `intervals` of the wrong length was silently recycled or gave an unrelated error.
* Loading chebpol no longer overwrites the `chebpol.threads` option when it is already set, for example in `.Rprofile`. The `CHEBPOL_THREADS` environment variable is only used when the option is unset.
* `evalongrid()` and `ipol()` no longer overwrite an argument vector that the function has kept. Previously one vector was reused for every grid point, so a function that kept its argument, for example to record the points it was called with, saw every copy change to the last grid point.
* `evalongrid()` and `ipol()` now call the function once at each grid point. Previously it was called twice at the first grid point, so a function that draws random numbers, counts its calls or has other side effects did so once more than the grid has points.
* `evalongrid()` now stops with an error when a grid dimension has no points. Previously it crashed R.
* `evalongrid()` and `ipol()` now stop with an error when `grid` has more points than R can store, e.g. four dimensions of 65,536 points each. Previously the number of points overflowed and R could crash.
* `evalongrid()` and `ipol()` now stop with an error when `grid` has no dimensions, e.g. `grid = list()`. Previously `ipol()` returned an interpolant without variables for several methods, and other methods failed with "'dim' cannot be of length 0".
* `evalongrid()` and `ipol()` now accept functions that return integer or logical values, and stop with an error when the length of the returned value changes between grid points.
* `ipol()` interpolants created with an explicit `blend` now keep their `arity` and `domain` attributes. Previously setting the default blend dropped them.
* `ipol()` now stops with an error when `grid` contains missing or infinite values, or values that are not numbers, such as strings. Previously these failed with "missing value where TRUE/FALSE needed" or gave an interpolant that returned `NaN`.
* `ipol()` with `method = "chebyshev"` or `"uniform"` now stops with an error when the function does not return a single number. Previously a function returning a vector gave an interpolant with one variable too many, treating the position in the returned vector as an extra dimension.
* `ipol()` now stops with an error when `val` is neither numeric values nor a function, such as a function name given as a string. Previously such input failed with unrelated errors, and a list of values was accepted by some methods.
* `ipol()` now stops with an error when `blend` is not supported by the chosen method. Previously `blend = "mean"` was accepted for `method = "multilinear"` and `"oldstalker"` and silently treated as linear blending.
* `ipol()` interpolants created with `method = "multilinear"` or `"oldstalker"` no longer return huge or `NaN` values just outside the grid with a blend other than `"linear"`, e.g. about 3.6e8 at -0.01 for `sin()` on `[0, 1]` with `blend = "sigmoid"`. Outside the grid these blends now give the value at the nearest edge for `"multilinear"`, and continue linearly from the edge for `"oldstalker"`. Values inside the grid, and with `blend = "linear"`, are unchanged.
* `ipol()` interpolants created with `method = "multilinear"` now return `NA` for `NA` or `NaN` input. Previously the evaluation never returned. `ipol()` also stops with an error for a multilinear grid dimension with a single point, which was read out of bounds.
* `ipol()` interpolants created with `method = "simplexlinear"` no longer return `NaN` with `blend = "square"` or `blend = "sigmoid"`. Square blending now uses the value at the vertex with the largest barycentric coordinate, and coordinates are clamped to `[0, 1]` before blending. `NA` or `NaN` input now gives `NA` for every blend.
* `ipol()` with `method = "fh"` now stops with an error when the degree `k` is negative, not a whole number or not smaller than the number of grid points in a dimension. Before, the interpolant returned `NaN` between the grid points, and `k = 1.5` was silently used as 1.
* `ipol()` with `method = "fh"` now stops with an error when `k` has neither length 1 nor one value per grid dimension. Previously a shorter `k` was silently repeated, or failed with "Length of grid (3) should equal length of d (2)".
* `ipol()` with `method = "fh"` now forwards extra arguments in `...` to the function being interpolated.
* `ipol()` with `method = "fh"` now stops with an error when the dimensions of the values do not match the grid lengths. Previously the interpolant read past the end of the grid and returned meaningless values.
* `ipol()` with `method = "fh"` no longer depends on the scale of the grid. Previously any point within 2e-15 of a grid point got the value at that grid point, which made the interpolant inaccurate on grids spaced less than about 1e-12 apart.
* `ipol()` with `method = "uniform"` now accepts a one-dimensional `grid` given as a bare vector, as the other grid methods do. Previously the interpolant returned `NaN`s.
* `ipol()` with `method = "general"` now handles a grid dimension with a single point, along which the interpolant is constant, as `method = "fh"` does. Previously the interpolant returned `NA` everywhere.
* `ipol()` with `method = "uniform"` now handles a reversed interval such as `c(1, 0)` when given a function. Previously the interpolant was silently wrong.
* `ipol()` with `method = "uniform"` now accepts a single interval given as a vector `c(a, b)` together with values. It already did when given a function.
* `ipol()` with `method = "simplexlinear"` now forwards extra arguments in `...` to the function being interpolated.
* `ipol()` with `method = "simplexlinear"` now stops with an error when the number of values differs from the number of knots. Previously the interpolant read past the end of the values and returned meaningless values or `NA`.
* `ipol()` with `method = "simplexlinear"` now accepts a function that returns integer or logical values, as the other methods do. Previously the interpolant failed with "REAL() can only be applied to a 'numeric', not a 'integer'" when it was evaluated.
* `ipol()` with `method = "simplexlinear"` now accepts one-dimensional knots given as a vector, as `method = "polyharmonic"` does. Previously it failed with "argument is of length zero".
* `ipol()` with `method = "polyharmonic"` accepts `normalize` as a vector `c(a, b)` or a two-column matrix again. Both forms failed with "the condition has length > 1" on R >= 4.2, and the knots were never transformed, so the fit did not match the evaluation.
* `ipol()` interpolants created with `method = "polyharmonic"` now report the range of the knots in their `domain` attribute. With normalization it was the range after scaling, `[0, 1]`.
* `ipol()` with `method = "polyharmonic"` no longer drops extra arguments in `...` for the function being interpolated when neither `normalize` nor `nowarn` is given.
* `ipol()` with `method = "polyharmonic"` now computes the size of the knots-by-knots distance matrix correctly for 46,341 or more knots. Previously the size computation overflowed, giving "negative length vectors are not allowed" or, from 65,536 knots, possibly writing outside the allocated memory. Fitting that many knots still needs enough memory for the matrix, about 17 GB at 46,341 knots.
* `ipol()` with `method = "polyharmonic"` now fits knots on a line or plane, such as a single knot or knots sharing a coordinate, and knots with a large offset or with coordinates on very different scales. Previously these returned `NaN` or fell back to an inaccurate least squares fit with a warning.
* `ipol()` with `method = "polyharmonic"` now fits a spline with a single knot. Previously it failed with "invalid 'nrow' value".
* `ipol()` with `method = "polyharmonic"` now stops with an error when `k` is zero, missing, infinite or not a single number. Previously `k = 0` silently fitted a spline with a logarithmic kernel that could be far off between the knots, and the other cases failed with "missing value where TRUE/FALSE needed" or "the condition has length > 1".
* `ipol()` with `method = "polyharmonic"` now stops with an error when `nowarn` is not `TRUE` or `FALSE`. Previously `nowarn = NA` failed with "missing value where TRUE/FALSE needed", and only when the fit fell back to least squares.
* `ipol()` with `method = "polyharmonic"` now stops with an error when the number of values differs from the number of knots. Previously it warned about falling back to a least squares fit and suggested `normalize = TRUE`.
* `ipol()` with `method = "stalker"`, `"hstalker"` or `"oldstalker"` now stops with an error when a grid dimension has a single point. Previously the spline read past the end of the values and the grid.
* `ipol()` interpolants created with `method = "hstalker"` no longer return `NaN` at a grid point with `blend = "mean"` when two neighboring values are equal.
* `ipol()` with `method = "hstalker"` no longer depends on the scale of the values or of the grid. Previously small differences between neighboring values were treated as constant or linear, and could produce a narrow spike between grid points.
* `ipol()` with `method = "stalker"` no longer depends on the scale of the values. Previously small differences between neighboring values were treated as constant or linear, so interpolating small values, e.g. `1e-6 * f`, gave a different curve than interpolating `f` and scaling the result.
* `ipol()` with `method = "stalker"` on a uniform grid now gives the same spline whether or not chebpol was built with GSL. Previously a build without GSL used degree 2 next to every local extremum, where a build with GSL lowers the degree to limit the overshoot. The warning that a non-uniform grid needs GSL no longer depends on the scale of the grid either: a non-uniform grid spaced less than about 1e-5 apart gave no warning, and a uniform grid with coordinates around 1e12 warned.
* `ipol()` with `method = "oldstalker"` and `k = NA` no longer depends on the scale of the grid. Previously a non-uniform grid spaced less than about 1e-8 apart was treated as uniform, and the interpolant could overshoot the data by orders of magnitude.
* `ipol()` with `method = "oldstalker"` and `k = 0` no longer depends on the scale of the values or of the grid. Previously values below about 1e-9 or a grid spaced less than about 1e-6 apart gave a different interpolant, which could overshoot the data by orders of magnitude at a local extremum.
* `ipol()` interpolants created with `method = "simplexlinear"` now stop with an error when `epol` is not `TRUE` or `FALSE`. Previously `epol = NA` extrapolated as if it were `TRUE`.
* `ipol()` interpolants created with `method = "stalker"`, `"hstalker"` or `"simplexlinear"` now stop with an error when given a matrix with the wrong number of rows. Previously they returned values silently, and could read past the end of the input.
* `ipol()` interpolants created with `method = "oldstalker"` now report both ranks when given input of the wrong dimension.

## Documentation

* The vignettes are now R Markdown documents rendered to HTML instead of Sweave documents rendered to PDF.

# chebpol 2.0

* The stalker spline now works in N dimensions.
* A new hyperbolic stalker has also been introduced.
