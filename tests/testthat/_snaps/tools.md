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

