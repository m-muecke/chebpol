# chebpol (development version)

## Breaking changes

* chebpol now requires R >= 4.1.0.
* `chebappx()`, `chebappxf()`, `chebappxg()`, `chebappxgf()`, `fhappx()`, `mlappx()`, `polyh()`, `slappx()`, `ucappx()` and `ucappxf()`, deprecated since version 1.7, have been removed. Use `ipol()` instead.

## Bug fixes

* chebpol compiles again with current versions of R. FFTW is now found when it is installed outside the default include path.
* `ipol()` with `method = "fh"` now forwards extra arguments in `...` to the function being interpolated.
* `ipol()` with `method = "polyharmonic"` no longer drops extra arguments in `...` for the function being interpolated when neither `normalize` nor `nowarn` is given.
* `ipol()` interpolants created with `method = "oldstalker"` now report both ranks when given input of the wrong dimension.

## Documentation

* The vignettes are now R Markdown documents rendered to HTML instead of Sweave documents rendered to PDF.

# chebpol 2.0

* The stalker spline now works in N dimensions.
* A new hyperbolic stalker has also been introduced.
