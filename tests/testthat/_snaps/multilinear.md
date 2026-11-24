# ipol() rejects a blend that multilinear interpolation does not support

    Code
      ipol(f, grid = grid, method = "multilinear", blend = "mean")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "linear", "cubic", "sigmoid", "parodic", "square"

