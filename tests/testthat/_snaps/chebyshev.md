# chebcoef() errors when dct is not TRUE or FALSE

    Code
      chebcoef(1:3, dct = NA)
    Condition
      Error in `chebcoef()`:
      ! `dct` must be `TRUE` or `FALSE`.
    Code
      chebcoef(1:3, dct = c(TRUE, FALSE))
    Condition
      Error in `chebcoef()`:
      ! `dct` must be `TRUE` or `FALSE`.

# ipol() errors on malformed Chebyshev intervals

    Code
      ipol(val, intervals = "a", method = "chebyshev")
    Condition
      Error in `chebappx.real()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.
    Code
      ipol(val, intervals = list(c(0, 1)), method = "chebyshev")
    Condition
      Error in `chebappx.real()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.
    Code
      ipol(val, intervals = list(c(0, 1), c(0, 1, 2)), method = "chebyshev")
    Condition
      Error in `chebappx.real()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.

