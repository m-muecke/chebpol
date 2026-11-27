# printing an interpolant summarizes it

    Code
      ipol(function(x) sum(x), grid = grid, method = "multilinear", blend = "cubic")
    Output
      <chebpol_interpolant> multilinear, 2 variables
      Domain: [0, 1] x [0, 2]
      Defaults: blend = "cubic"
    Code
      ipol(function(x) sum(x), dims = 6, method = "chebyshev")
    Output
      <chebpol_interpolant> chebyshev, 1 variable
      Domain: [-1, 1]
    Code
      ipol(function(x) sum(x), grid = grid, method = "oldstalker", k = 1.5)
    Output
      <chebpol_interpolant> oldstalker, 2 variables
      Domain: [0, 1] x [0, 2]
      Defaults: degree = c(1.5, 1.5), blend = "linear"

# printing a simplex linear interpolant shows the range of the knots

    Code
      ipol(function(x) sum(x), knots = knots, method = "simplexlinear")
    Output
      <chebpol_interpolant> simplexlinear, 2 variables
      Domain: [0, 2] x [0, 1]
      Defaults: epol = FALSE, blend = "linear"

