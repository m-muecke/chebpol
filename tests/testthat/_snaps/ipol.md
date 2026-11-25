# ipol() errors on a grid without dimensions

    Code
      ipol(1, grid = list(), method = "multilinear")
    Condition
      Error in `ipol()`:
      ! `grid` must have at least one dimension.
    Code
      ipol(sum, grid = list(), method = "fh")
    Condition
      Error in `ipol()`:
      ! `grid` must have at least one dimension.

