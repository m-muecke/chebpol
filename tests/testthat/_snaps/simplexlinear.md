# simplex linear interpolants error when epol is not TRUE or FALSE

    Code
      sl(c(2, 2), epol = NA)
    Condition
      Error in `sl()`:
      ! `epol` must be `TRUE` or `FALSE`.
    Code
      sl(c(2, 2), epol = "yes")
    Condition
      Error in `sl()`:
      ! `epol` must be `TRUE` or `FALSE`.

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
      Error in `sl()`:
      ! `x` must have 2 rows, one per dimension, not 1.
    Code
      sl(matrix(0, 3, 2))
    Condition
      Error in `sl()`:
      ! `x` must have 2 rows, one per dimension, not 3.

# ipol() errors when values and knots differ in number for simplex linear interpolation

    Code
      ipol(1:3, knots = knots, method = "simplexlinear")
    Condition
      Error in `ipol_simplexlinear()`:
      ! `val` must have one value per knot (10), not 3.

# ipol() errors for simplex linear interpolation without geometry

    Code
      ipol(sum, knots = matrix(runif(20), 2L), method = "simplexlinear")
    Condition
      Error in `ipol_simplexlinear()`:
      ! method = "simplexlinear" requires the geometry package.

