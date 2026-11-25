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

# evalongrid() errors when the function value is not usable

    Code
      evalongrid(function(x) "a", 3)
    Condition
      Error in `evalongrid()`:
      ! `fun` must return a numeric value.
    Code
      evalongrid(growing, 4)
    Condition
      Error in `evalongrid()`:
      ! `fun` must return a value of length 1 at every grid point, not 2.

# evalongrid() errors on an empty grid dimension

    Code
      evalongrid(sum, grid = list(numeric(0), 1:3))
    Condition
      Error in `evalongrid()`:
      ! `grid` dimension 1 must have at least one point.
    Code
      evalongrid(sum, dims = c(2, 0))
    Condition
      Error in `evalongrid()`:
      ! `grid` dimension 2 must have at least one point.
    Code
      evalongrid(sum, grid = list())
    Condition
      Error in `evalongrid()`:
      ! `grid` must have at least one dimension.
    Code
      evalongrid(sum, dims = integer(0))
    Condition
      Error in `evalongrid()`:
      ! `grid` must have at least one dimension.

# evalongrid() errors on a grid with too many points

    Code
      evalongrid(sum, grid = list(axis, axis, axis, axis))
    Condition
      Error in `evalongrid()`:
      ! `grid` has too many points.
    Code
      evalongrid(function(x) numeric(8), grid = rep(list(axis[1:1024]), 5))
    Condition
      Error in `evalongrid()`:
      ! `fun` returns too many values on `grid` to store them.

