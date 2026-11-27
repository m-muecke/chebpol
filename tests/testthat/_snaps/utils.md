# interpolants accept a point, a matrix, or a vector of several points

    Code
      out <- ml(c(m))
    Condition
      Warning in `ml()`:
      Treating `x` of length 6 as a 2 x 3 matrix, one point per column.

# interpolants error when the input length does not fit the number of variables

    Code
      ml(c(0.1, 0.2))
    Condition
      Error in `ml()`:
      ! `x` must have length 3 or be a matrix with 3 rows, not length 2.
    Code
      gen(c(0.1, 0.2))
    Condition
      Error in `gen()`:
      ! `x` must have length 3 or be a matrix with 3 rows, not length 2.
    Code
      gen(matrix(0, 2, 4))
    Condition
      Error in `gen()`:
      ! `x` must have 3 rows, one per dimension, not 2.
    Code
      uni(matrix(0, 2, 4))
    Condition
      Error in `uni()`:
      ! `x` must have 3 rows, one per dimension, not 2.

