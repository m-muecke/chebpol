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

# ipol() errors when values and knots differ in number

    Code
      ipol(1:3, knots = knots, method = "polyharmonic")
    Condition
      Error:
      ! `val` must have one value per knot (10), not 3.

