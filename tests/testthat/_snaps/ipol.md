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

# ipol() errors when val is neither values nor a function

    Code
      ipol("sum", grid = grid, method = "fh")
    Condition
      Error in `ipol()`:
      ! `val` must be numeric values or a function, not a character.
    Code
      ipol(as.list(1:6), grid = grid, method = "multilinear")
    Condition
      Error in `ipol()`:
      ! `val` must be numeric values or a function, not a list.
    Code
      ipol(data.frame(a = 1:6), dims = 6, method = "chebyshev")
    Condition
      Error in `ipol()`:
      ! `val` must be numeric values or a function, not a data.frame.

