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

# ipol() checks the grid of the grid methods

    Code
      ipol(sum, method = "fh")
    Condition
      Error in `ipol()`:
      ! `grid` must be given for method = "fh".
    Code
      ipol(sum, grid = list(1), method = "multilinear")
    Condition
      Error in `ipol()`:
      ! Each `grid` dimension must have at least two points for method = "multilinear".
    Code
      ipol(sum, grid = c(0, 2, 1), method = "general")
    Condition
      Error in `ipol()`:
      ! Each `grid` dimension must be strictly increasing or strictly decreasing.

