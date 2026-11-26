# ipol() errors on a malformed normalize

    Code
      ipol(f, knots = knots, method = "polyharmonic", normalize = c(1, 2, 3))
    Condition
      Error:
      ! `normalize` must be a logical, a vector `c(a, b)`, or a two-column matrix.
    Code
      ipol(f, knots = knots, method = "polyharmonic", normalize = matrix(1, 3L, 2L))
    Condition
      Error:
      ! The number of rows of `normalize` (3) must divide the number of dimensions (2).

# ipol() errors on an invalid polyharmonic k

    Code
      ipol(f, knots = knots, k = 0, method = "polyharmonic")
    Condition
      Error:
      ! `k` must be a positive whole number or a negative number.
    Code
      ipol(f, knots = knots, k = NA, method = "polyharmonic")
    Condition
      Error:
      ! `k` must be a positive whole number or a negative number.
    Code
      ipol(f, knots = knots, k = c(2, 3), method = "polyharmonic")
    Condition
      Error:
      ! `k` must be a positive whole number or a negative number.
    Code
      ipol(f, knots = knots, k = 2.5, method = "polyharmonic")
    Condition
      Error:
      ! A positive `k` must be a whole number, not 2.5.

# ipol() errors when values and knots differ in number

    Code
      ipol(1:3, knots = knots, method = "polyharmonic")
    Condition
      Error:
      ! `val` must have one value per knot (10), not 3.

# polyharmonic interpolation warns for the ipol() call when it falls back to least squares

    Code
      ph <- ipol(c(1, 2, 3, 4), knots = knots, method = "polyharmonic")
    Condition
      Warning in `ipol()`:
      Can't fit the spline exactly, falling back to a least squares fit. Try `normalize = TRUE`.

