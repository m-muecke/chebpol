ipol_polyharmonic <- function(val, knots, k = 2, normalize = NA, nowarn = FALSE, ...) {
  # Linear polyharmonic splines. Centres are columns in matrix knots. Function values in val.
  # k = 2 yields thin-plate splines.
  # k < 0 yields Gaussian kernel splines, with sigma^2 = -1 / k
  # we compute r^2, so powers etc. are adjusted for that case
  if (!is_flag(nowarn)) {
    stop("`nowarn` must be `TRUE` or `FALSE`.")
  }
  if (is.null(dim(knots))) {
    dim(knots) <- c(1L, length(knots))
  }
  if (is.function(val)) {
    val <- apply(knots, 2L, val, ...)
  }
  n_knots <- ncol(knots)
  n_dim <- nrow(knots)
  if (length(val) != n_knots) {
    stop(sprintf("`val` must have one value per knot (%d), not %d.", n_knots, length(val)))
  }
  domain <- knots_domain(knots)
  spans <- vapply(domain, diff, numeric(1L), USE.NAMES = FALSE)
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k == 0) {
    stop("`k` must be a positive whole number or a negative number.")
  }
  if (k > 0) {
    if (abs(round(k) - k) > sqrt(.Machine$double.eps)) {
      stop(sprintf("A positive `k` must be a whole number, not %s.", format(k)))
    }
    k <- as.integer(round(k))
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
    if (n_dim %% nrow(normalize) != 0L) {
      stop(sprintf(
        "The number of rows of `normalize` (%d) must divide the number of dimensions (%d).",
        nrow(normalize),
        n_dim
      ))
    }
    scale <- normalize[, 1L]
    offset <- normalize[, 2L]
    normfun <- crate(function(x) scale * (x - offset), scale, offset)
    knots <- normfun(knots)
    normalize <- TRUE
  } else if (normalize) {
    scale <- 1 / ifelse(spans > 0, spans, 1)
    offset <- apply(knots, 1L, min)
    normfun <- crate(function(x) scale * (x - offset), scale, offset)
    knots <- normfun(knots)
    normalize <- TRUE
  } else {
    normfun <- identity
    normalize <- FALSE
  }

  threads <- getOption("chebpol.threads")
  # The distance matrix is unreferenced, so c_rbf_kernel transforms it in place
  A <- .Call(c_rbf_kernel, .Call(c_sq_dists, knots, knots, threads), k, threads)
  # fit the linear term in the directions the centered knots span, so that knots on a line or plane and large
  # offsets give a regular, well-conditioned system
  mu <- rowMeans(knots)
  centered <- knots - mu
  s <- svd(centered)
  # a direction is flat when the knots spread less along it than 1e-10 of their largest spread or than the rounding
  # of the coordinates it combines, so that a coordinate on a much smaller scale than the others still counts
  size <- sqrt(n_knots) * drop(crossprod(abs(s$u), apply(abs(knots), 1L, max)))
  flat <- pmax(1e-10 * pmin(s$d[1L], size), max(dim(knots)) * .Machine$double.eps * size)
  spread <- s$d > flat
  U <- s$u[, spread, drop = FALSE]
  # orthonormal rows on the scale of the distance matrix keep the system well scaled whatever the spread
  balance <- norm(A, "F")
  if (balance == 0) {
    balance <- 1
  }
  B <- balance * rbind(1 / sqrt(n_knots), t(s$v[, spread, drop = FALSE]))
  L <- nrow(B)
  n <- n_knots + L
  mat <- matrix(0, n, n)
  mat[seq_len(n_knots), seq_len(n_knots)] <- A
  mat[n_knots + seq_len(L), seq_len(n_knots)] <- B
  mat[seq_len(n_knots), n_knots + seq_len(L)] <- t(B)
  # free the distance matrix before solve() makes its own copy of mat
  rm(A, B)
  rhs <- c(val, rep(0, L))
  coefs <- try(solve(mat, rhs), silent = TRUE)
  if (inherits(coefs, "try-error")) {
    if (!nowarn) {
      msg <- paste0(
        "Can't fit the spline exactly, falling back to a least squares fit.",
        if (normalize) "" else " Try `normalize = TRUE`."
      )
      warning(warningCondition(msg, call = sys.call(sys.parent())))
    }
    coefs <- stats::lm.fit(mat, rhs)$coefficients
    coefs[is.na(coefs)] <- 0
  }
  w <- coefs[seq_len(n_knots)]
  basis_coef <- balance * coefs[n_knots + seq_len(L)]
  lin <- drop(U %*% (basis_coef[-1L] / s$d[spread]))
  v <- c(basis_coef[1L] / sqrt(n_knots) - sum(lin * mu), lin)
  fun <- function(x, threads = getOption("chebpol.threads")) {
    x <- as_points(x, n_dim)
    .Call(c_eval_polyharmonic, normfun(x), knots, w, v, k, threads)
  }
  new_interpolant(crate(fun, normfun, knots, w, v, k, n_dim), method = "polyharmonic", arity = n_dim, domain = domain)
}
