# ipol() errors without knots for simplex linear interpolation

    Code
      ipol(f, method = "simplexlinear")
    Condition
      Error in `ipol()`:
      ! `knots` or `grid` must be given for method = "simplexlinear".

# simplex linear interpolants error on input of the wrong dimension

    Code
      sl(matrix(0, 1, 5))
    Condition
      Error:
      ! `x` must have 2 rows, not 1.
    Code
      sl(matrix(0, 3, 2))
    Condition
      Error:
      ! `x` must have 2 rows, not 3.

# ipol() errors when values and knots differ in number for simplex linear interpolation

    Code
      ipol(1:3, knots = knots, method = "simplexlinear")
    Condition
      Error in `slappx.real()`:
      ! `val` must have one value per knot (10), not 3.

# ipol() errors for simplex linear interpolation without geometry

    Code
      ipol(sum, knots = matrix(runif(20), 2L), method = "simplexlinear")
    Condition
      Error in `slappx.real()`:
      ! method = "simplexlinear" requires the geometry package.

