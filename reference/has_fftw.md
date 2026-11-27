# Check which optional libraries chebpol was built with

chebpol can be compiled with or without two optional C libraries, each
of which enables part of the package:

## Usage

``` r
has_fftw()

has_gsl()
```

## Value

`TRUE` if chebpol was built with the library, otherwise `FALSE`.

## Details

- FFTW computes Chebyshev coefficients with a fast DCT. Without it,
  [`chebcoef()`](https://m-muecke.github.io/chebpol/reference/chebcoef.md)
  and [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md)
  with `method = "chebyshev"` use a matrix product, which is not
  feasible for high-degree polynomials: a 100x100x100 approximation is
  possible, but not a one-dimensional 1000000.

- GSL solves for the degree of the stalker splines on a non-uniform
  grid. Without it,
  [`ipol()`](https://m-muecke.github.io/chebpol/reference/ipol.md) with
  `method = "stalker"` uses degree 2 there and warns, and
  `method = "oldstalker"` with `k = NA` stops with an error.

When chebpol is installed from source on Linux or macOS, FFTW and GSL
are used if pkg-config finds them. To build without one, set the
environment variable `CHEBPOL_DISABLE_FFTW` or `CHEBPOL_DISABLE_GSL` to
a non-empty value before installing, or pass `--disable-fftw` or
`--disable-gsl` in `configure.args`.

## Examples

``` r
has_fftw()
#> [1] TRUE
has_gsl()
#> [1] TRUE
```
