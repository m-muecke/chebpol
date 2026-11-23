# oldstalker interpolant errors on input of the wrong dimension

    Code
      st(matrix(0, 3, 4))
    Condition
      Error:
      ! Rank of input (3) does not match rank of spline (2)

# ipol() errors on invalid stalker input

    Code
      ipol(f, method = "stalker")
    Condition
      Error in `ipol()`:
      ! grid must be specified for stalker interpolation
    Code
      ipol(f, grid = list(c(1, 0, 2)), method = "stalker")
    Condition
      Error in `ipol()`:
      ! grid must be distinct ordered values
    Code
      ipol(f, grid = grid, method = "stalker", blend = "nope")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "linear", "cubic", "sigmoid", "parodic", "square", "mean"

