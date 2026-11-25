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

# ipol() errors when the values do not match the Floater-Hormann grid

    Code
      ipol(array(0, c(3, 4)), grid = grid, method = "fh")
    Condition
      Error in `fhappx.real()`:
      ! dimensions of values (3 x 4) must match the grid lengths (4 x 3)

