# chebknots() errors on malformed dims and intervals

    Code
      chebknots(2.5)
    Condition
      Error in `chebknots()`:
      ! `dims` must be non-negative whole numbers.
    Code
      chebknots(c(2, 2), list(c(0, 1)))
    Condition
      Error in `chebknots()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `dims`.
    Code
      chebknots(2, list(c(0, 1), c(5, 6)))
    Condition
      Error in `chebknots()`:
      ! `intervals` must be a list of 1 `c(lower, upper)` pairs, one per dimension of `dims`.
    Code
      chebknots(3, c(0, 1, 2))
    Condition
      Error in `chebknots()`:
      ! `intervals` must be a list of 1 `c(lower, upper)` pairs, one per dimension of `dims`.

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
      Error in `ipol_chebyshev()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.
    Code
      ipol(val, intervals = list(c(0, 1)), method = "chebyshev")
    Condition
      Error in `ipol_chebyshev()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.
    Code
      ipol(val, intervals = list(c(0, 1), c(0, 1, 2)), method = "chebyshev")
    Condition
      Error in `ipol_chebyshev()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `val`.
    Code
      ipol(sum, dims = 3, intervals = list(c(0, 1), c(0, 2)), method = "chebyshev")
    Condition
      Error in `chebknots()`:
      ! `intervals` must be a list of 1 `c(lower, upper)` pairs, one per dimension of `dims`.
    Code
      ipol(sin, dims = 2.5, method = "chebyshev")
    Condition
      Error in `chebknots()`:
      ! `dims` must be non-negative whole numbers.

# chebeval() errors on malformed intervals

    Code
      chebeval(c(0.1, 0.2), coef, intervals = list(c(0, 1)))
    Condition
      Error in `chebeval()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `coef`.
    Code
      chebeval(c(0.1, 0.2), coef, intervals = list(c(0, 1), c(0, 1, 2)))
    Condition
      Error in `chebeval()`:
      ! `intervals` must be a list of 2 `c(lower, upper)` pairs, one per dimension of `coef`.

# chebeval() errors on empty coefficients

    Code
      chebeval(0, numeric())
    Condition
      Error in `chebeval()`:
      ! `coef` must have at least one element.
    Code
      chebeval(c(0, 0), matrix(numeric(), nrow = 0, ncol = 2))
    Condition
      Error in `chebeval()`:
      ! `coef` must have at least one element.

