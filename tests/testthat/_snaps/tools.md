# evalongrid() errors when the function value is not usable

    Code
      evalongrid(function(x) "a", 3)
    Condition
      Error in `evalongrid()`:
      ! fun must return a numeric value
    Code
      evalongrid(growing, 4)
    Condition
      Error in `evalongrid()`:
      ! fun must return a value of length 1, not 2

