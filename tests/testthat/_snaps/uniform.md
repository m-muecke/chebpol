# ipol() errors on malformed uniform intervals

    Code
      ipol(f, dims = c(5, 5), intervals = list(c(0, 1)), method = "uniform")
    Condition
      Error in `ipol_uniform()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.
    Code
      ipol(f, dims = c(5, 5), intervals = list(c(0, 1), c(0, 2), c(0, 3)), method = "uniform")
    Condition
      Error in `ipol_uniform()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.
    Code
      ipol(array(0, c(3, 4)), intervals = list(c(0, 1)), method = "uniform")
    Condition
      Error in `ipol_uniform()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.
    Code
      ipol(array(0, c(3, 4)), intervals = "a", method = "uniform")
    Condition
      Error in `ipol_uniform()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.

