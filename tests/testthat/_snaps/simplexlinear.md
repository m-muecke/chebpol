# ipol() errors without knots for simplex linear interpolation

    Code
      ipol(f, method = "simplexlinear")
    Condition
      Error in `ipol()`:
      ! knots must be specified for simplex linear interpolation

# simplex linear interpolants error on input of the wrong dimension

    Code
      sl(matrix(0, 1, 5))
    Condition
      Error:
      ! Dimension of input (1) must match dimension of interpolant (2)
    Code
      sl(matrix(0, 3, 2))
    Condition
      Error:
      ! Dimension of input (3) must match dimension of interpolant (2)

