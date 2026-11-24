# ipol() rejects a Floater-Hormann degree outside the valid range

    Code
      ipol(f, grid = grid, method = "fh", k = -1)
    Condition
      Error in `ipol()`:
      ! k must be between 0 and one less than the number of grid points in each dimension
    Code
      ipol(f, grid = grid, method = "fh", k = 4)
    Condition
      Error in `ipol()`:
      ! k must be between 0 and one less than the number of grid points in each dimension
    Code
      ipol(f, grid = grid, method = "fh", k = c(6, 3))
    Condition
      Error in `ipol()`:
      ! k must be between 0 and one less than the number of grid points in each dimension

