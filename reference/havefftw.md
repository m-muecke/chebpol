# Check whether chebpol uses FFTW

It is possible to compile chebpol without FFTW. If this is done, it will
not be feasible to use high-degree Chebyshev polynomials. I.e. a
100x100x100 approximation will be possible, but not a one-dimensional
1000000. This function checks whether chebpol uses FFTW.

## Usage

``` r
havefftw()
```

## Value

Returns `TRUE` if chebpol uses FFTW. Otherwise `FALSE`.
