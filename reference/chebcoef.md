# Compute Chebyshev-coefficients given values on a Chebyshev grid

Compute the multivariate Chebyshev-coefficients, given values on a
Chebyshev grid.

## Usage

``` r
chebcoef(val, dct = FALSE)
```

## Arguments

- val:

  An `array` of function values on a Chebyshev grid. The `dim`-attribute
  must be appropriately set. If not set, it is assumed to be
  one-dimensional.

- dct:

  Logical. Since the Chebyshev coefficients are closely related to the
  DCT-II transform of `val`, the non-normalized real-even DCT-II
  coefficients may be retrieved instead. I.e. those from `FFTW_REDFT10`
  in each dimension. This is not used anywhere in the package, it is
  merely provided as a convenience for those who might need it.

## Value

An array of Chebyshev-coefficients for an interpolating
Chebyshev-polynomial.

## Details

If `val` has no `dim`-attribute, it is assumed to be one-dimensional
with length `length(val)`.

If chebpol was compiled without FFTW, running `chebcoef()` on large
grids may be slow and memory-demanding.

## See also

[`havefftw()`](https://m-muecke.github.io/chebpol/reference/havefftw.md)

## Examples

``` r
# Coefficients for a 2x3x4 grid
a <- array(rnorm(24), dim = c(2, 3, 4))
chebcoef(a)
#> , , 1
#> 
#>            [,1]       [,2]       [,3]
#> [1,] -0.2473444 -0.2781946  0.1888079
#> [2,] -0.1302599  0.6271207 -0.1457199
#> 
#> , , 2
#> 
#>            [,1]       [,2]      [,3]
#> [1,] -0.1645531 -0.5889327 0.5128264
#> [2,] -0.8253577 -1.2287445 0.3194804
#> 
#> , , 3
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  0.1230125 -0.1615874 0.07603228
#> [2,] -0.6647776 -0.2024659 0.52033074
#> 
#> , , 4
#> 
#>            [,1]      [,2]      [,3]
#> [1,] 0.02470929 0.2486606 0.5281751
#> [2,] 0.71952494 0.4979239 0.2994270
#> 
```
