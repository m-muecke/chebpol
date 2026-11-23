\name{NEWS}
\title{chebpol news}
\section{Changes in version 2.1-2.9000}{
  \itemize{
    \item The functions \code{chebappx()}, \code{chebappxf()}, \code{chebappxg()}, \code{chebappxgf()}, \code{fhappx()}, \code{mlappx()}, \code{polyh()}, \code{slappx()}, \code{ucappx()} and \code{ucappxf()}, deprecated since version 1.7, have been removed. Use \code{ipol()} instead.
    \item \code{ipol()} with \code{method = "fh"} now forwards extra arguments in \code{...} to the function being interpolated.
  }
}
\section{Changes in version 2.0}{
  \itemize{
    \item The stalker spline now works in N dimensions.
    \item A new hyperbolic stalker has also been introduced.
  }
}
