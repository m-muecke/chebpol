# ipol() rejects a blend that multilinear interpolation does not support

    Code
      ipol(f, grid = grid, method = "multilinear", blend = "mean")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "linear", "cubic", "sigmoid", "parodic", "square"

# ipol() rejects a multilinear grid dimension with a single point

    Code
      ipol(c(1, 2, 3), grid = list(5, c(1, 2, 3)), method = "multilinear")
    Condition
      Error in `ipol()`:
      ! each grid dimension must have at least two points for multi linear interpolation

