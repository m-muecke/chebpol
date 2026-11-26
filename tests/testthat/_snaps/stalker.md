# oldstalker interpolant errors on input of the wrong dimension

    Code
      st(matrix(0, 3, 4))
    Condition
      Error in `st()`:
      ! `x` must have 2 rows, one per dimension, not 3.

# ipol() errors on invalid stalker input

    Code
      ipol(f, method = "stalker")
    Condition
      Error in `ipol()`:
      ! `grid` must be given for method = "stalker".
    Code
      ipol(f, grid = list(c(1, 0, 2)), method = "stalker")
    Condition
      Error in `ipol()`:
      ! Each `grid` dimension must be strictly increasing or strictly decreasing.
    Code
      ipol(f, grid = grid, method = "stalker", blend = "nope")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "cubic", "linear", "sigmoid", "parodic", "square", "mean"

# ipol() errors when the values do not match the stalker grid

    Code
      ipol(1:5, grid = grid, method = "stalker")
    Condition
      Error in `ipol_stalker()`:
      ! `val` must have one value per grid point (9), not 5.
    Code
      ipol(1:5, grid = grid, method = "hstalker")
    Condition
      Error in `ipol_hstalker()`:
      ! `val` must have one value per grid point (9), not 5.
    Code
      ipol(1:5, grid = grid, method = "oldstalker")
    Condition
      Error in `ipol_oldstalker()`:
      ! `val` must have one value per grid point (9), not 5.

# ipol() errors on a stalker grid dimension with a single point

    Code
      ipol(sum, grid = grid, method = "stalker")
    Condition
      Error in `ipol()`:
      ! Each `grid` dimension must have at least two points for method = "stalker".
    Code
      ipol(sum, grid = grid, method = "hstalker")
    Condition
      Error in `ipol()`:
      ! Each `grid` dimension must have at least two points for method = "hstalker".
    Code
      ipol(sum, grid = grid, method = "oldstalker")
    Condition
      Error in `ipol()`:
      ! Each `grid` dimension must have at least two points for method = "oldstalker".

# stalker interpolants error on input of the wrong dimension

    Code
      st(matrix(0, 1, 5))
    Condition
      Error in `st()`:
      ! `x` must have 2 rows, one per dimension, not 1.
    Code
      st(matrix(0, 3, 2))
    Condition
      Error in `st()`:
      ! `x` must have 2 rows, one per dimension, not 3.
    Code
      hst(matrix(0, 1, 5))
    Condition
      Error in `hst()`:
      ! `x` must have 2 rows, one per dimension, not 1.
    Code
      hst(matrix(0, 3, 2))
    Condition
      Error in `hst()`:
      ! `x` must have 2 rows, one per dimension, not 3.

