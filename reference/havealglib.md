# Check whether chebpol has the ALGLIB library

If ALGLIB was available at compile time, it can be used for compact
support radial basis function interpolations on scattered data. This
function checks whether ALGLIB is available.

## Usage

``` r
havealglib()
```

## Value

Returns TRUE if ALGLIB is available. Otherwise FALSE.
