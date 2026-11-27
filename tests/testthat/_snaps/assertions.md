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

