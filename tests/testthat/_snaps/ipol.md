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

# ipol() errors on a grid with missing values

    Code
      ipol(1:3, grid = c(0, NA, 1), method = "multilinear")
    Condition
      Error in `ipol()`:
      ! `grid` must not contain missing values.
    Code
      ipol(sum, grid = list(0:2, c(0, NaN, 1)), method = "fh")
    Condition
      Error in `ipol()`:
      ! `grid` must not contain missing values.

# ipol() errors when the function does not return a single number

    Code
      ipol(function(x) c(x, 1), dims = 3, method = "chebyshev")
    Condition
      Error in `ipol_chebyshev_fun()`:
      ! The function `val` must return a single number, not 2 values.
    Code
      ipol(function(x) c(x, 1), dims = c(3, 4), method = "uniform")
    Condition
      Error in `ipol_uniform_fun()`:
      ! The function `val` must return a single number, not 3 values.
    Code
      ipol(function(x) numeric(0), dims = 3, method = "chebyshev")
    Condition
      Error in `ipol_chebyshev_fun()`:
      ! The function `val` must return a single number, not 0 values.

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

