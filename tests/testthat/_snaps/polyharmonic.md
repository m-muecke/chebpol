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

# ipol() errors when values and knots differ in number

    Code
      ipol(1:3, knots = knots, method = "polyharmonic")
    Condition
      Error:
      ! length of values (3) must equal the number of knots (10)

