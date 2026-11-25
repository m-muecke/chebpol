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

