polyh.real <- function(val, knots, k = 2, normalize = NA, nowarn = FALSE, ...) {
  # Linear polyharmonic splines. Centres are columns in matrix knots. Function values in val.
  # k = 2 yields thin-plate splines.
  # k < 0 yields Gaussian kernel splines, with sigma^2 = -1 / k
  # we compute r^2, so powers etc. are adjusted for that case
  if (is.null(dim(knots))) {
    dim(knots) <- c(1L, length(knots))
  }
  if (is.function(val)) {
    val <- apply(knots, 2L, val, ...)
  }
  N <- ncol(knots)
  M <- nrow(knots)
  if (length(val) != N) {
    stop(sprintf("`val` must have one value per knot (%d), not %d.", N, length(val)))
  }
  domain <- knots_domain(knots)
  spans <- vapply(domain, diff, numeric(1L), USE.NAMES = FALSE)
  if (k > 0) {
    if (abs(round(k) - k) > sqrt(.Machine$double.eps)) {
      stop(sprintf("A positive `k` must be a whole number, not %s.", format(k)))
    } else {
      k <- as.integer(round(k))
    }
  }

  if (length(normalize) == 1L && is.na(normalize)) {
    normalize <- min(knots) < 0 || max(knots) > 1
  } else if (!is.matrix(normalize) && length(normalize) == 2L) {
    dim(normalize) <- c(1L, 2L)
  }
  if (!is.matrix(normalize) && length(normalize) != 1L) {
    stop("`normalize` must be a logical, a vector `c(a, b)`, or a two-column matrix.")
  }
  if (is.matrix(normalize)) {
    if (ncol(normalize) != 2L) {
      stop("`normalize` must be a logical, a vector `c(a, b)`, or a two-column matrix.")
    }
    if (M %% nrow(normalize) != 0L) {
      stop(sprintf(
        "The number of rows of `normalize` (%d) must divide the number of dimensions (%d).",
        nrow(normalize),
        M
      ))
    }
    wa <- normalize[, 1L]
    wb <- normalize[, 2L]
    normfun <- crate(function(x) wa * (x - wb), wa, wb)
    knots <- normfun(knots)
    normalize <- TRUE
  } else if (normalize) {
    wa <- 1 / ifelse(spans > 0, spans, 1)
    wb <- apply(knots, 1L, min)
    normfun <- crate(function(x) wa * (x - wb), wa, wb)
    knots <- normfun(knots)
    normalize <- TRUE
  } else {
    normfun <- identity
    normalize <- FALSE
  }

  threads <- getOption("chebpol.threads")
  # The distance matrix is unreferenced, so C_phifunc transforms it in place
  A <- .Call(C_phifunc, .Call(C_sqdiffs, knots, knots, threads), k, threads)
  # drop linearly dependent rows of the linear term, e.g. when the knots lie on a line or a plane
  B <- rbind(1, knots)
  # the default tolerance treats a small spread around a large offset as constant
  q <- qr(t(B), tol = 1e-10)
  keep <- sort(q$pivot[seq_len(q$rank)])
  B <- B[keep, , drop = FALSE]
  L <- nrow(B)
  n <- N + L
  mat <- matrix(0, n, n)
  mat[seq_len(N), seq_len(N)] <- A
  mat[N + seq_len(L), seq_len(N)] <- B
  mat[seq_len(N), N + seq_len(L)] <- t(B)
  # free the distance matrix before solve() makes its own copy of mat
  rm(A, B)
  rhs <- c(val, rep(0, L))
  wv <- try(solve(mat, rhs), silent = TRUE)
  if (inherits(wv, "try-error")) {
    if (!nowarn) {
      msg <- paste0(
        "Can't fit the spline exactly, falling back to a least squares fit.",
        if (!normalize) " Try `normalize = TRUE`." else ""
      )
      warning(warningCondition(msg, call = sys.call(sys.parent())))
    }
    wv <- stats::lm.fit(mat, rhs)$coefficients
    wv[is.na(wv)] <- 0
  }
  w <- wv[seq_len(N)]
  v <- numeric(M + 1L)
  v[keep] <- wv[N + seq_len(L)]
  fun <- function(x, threads = getOption("chebpol.threads")) {
    x <- as_points(x, M)
    .Call(C_evalpolyh, normfun(x), knots, w, v, k, threads, NULL)
  }
  new_interpolant(crate(fun, normfun, knots, w, v, k, M), method = "polyharmonic", arity = M, domain = domain)
}
