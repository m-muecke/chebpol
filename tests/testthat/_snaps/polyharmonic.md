# ipol() errors on a malformed normalize

    Code
      ipol(f, knots = knots, method = "polyharmonic", normalize = c(1, 2, 3))
    Condition
      Error:
      ! normalize must be a logical, a vector c(a, b), or a matrix with two columns
    Code
      ipol(f, knots = knots, method = "polyharmonic", normalize = matrix(1, 3L, 2L))
    Condition
      Error:
      ! normalize must be a n x 2 matrix with 2 divisible by n (3)

