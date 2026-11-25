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
      ! 'arg' should be one of "cubic", "linear", "sigmoid", "parodic", "square", "mean"

# ipol() errors on a stalker grid dimension with a single point

    Code
      ipol(sum, grid = grid, method = "stalker")
    Condition
      Error in `ipol()`:
      ! each grid dimension must have at least two points for stalker interpolation
    Code
      ipol(sum, grid = grid, method = "hstalker")
    Condition
      Error in `ipol()`:
      ! each grid dimension must have at least two points for stalker interpolation
    Code
      ipol(sum, grid = grid, method = "oldstalker")
    Condition
      Error in `ipol()`:
      ! each grid dimension must have at least two points for stalker interpolation

# stalker interpolants error on input of the wrong dimension

    Code
      st(matrix(0, 1, 5))
    Condition
      Error:
      ! Rank of input (1) does not match rank of spline (2)
    Code
      st(matrix(0, 3, 2))
    Condition
      Error:
      ! Rank of input (3) does not match rank of spline (2)
    Code
      hst(matrix(0, 1, 5))
    Condition
      Error:
      ! Rank of input (1) does not match rank of spline (2)
    Code
      hst(matrix(0, 3, 2))
    Condition
      Error:
      ! Rank of input (3) does not match rank of spline (2)

