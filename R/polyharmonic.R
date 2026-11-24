polyh.real <- function(val, knots, k = 2, normalize = NA, nowarn = FALSE, ...) {
  # Linear polyharmonic splines. Centres are columns in matrix knots. Function values in val.
  # Quite slow for ncol(knots) > 3000 or so.  k=2 yields thin-plate splines.
  # There exist faster evaluation methods for dimensions <= 4
  # k < 0 yields Gaussian kernel splines, with sigma^2 = -1/k
  # we compute r^2, so powers etc. are adjusted for that case
  # Hmm, I should take a look at http://dx.doi.org/10.1016/j.jat.2012.11.008 for the unit ball
  # perhaps some normalization should be added?
  # I also need to look at fast summation with NFFT
  # https://www-user.tu-chemnitz.de/~potts/nfft/fastsum.php
  # Likewise its interpolation applications which I don't understand yet.
  if (is.null(dim(knots))) {
    dim(knots) <- c(1L, length(knots))
  }
  if (is.function(val)) {
    val <- apply(knots, 2L, val, ...)
  }
  N <- ncol(knots)
  M <- nrow(knots)
  if (k > 0) {
    if (abs(round(k) - k) > sqrt(.Machine$double.eps)) {
      stop(sprintf("k is positive, but not an integer: %.17f", k))
    } else {
      k <- as.integer(round(k))
    }
  }

  if (is.na(normalize)) {
    normalize <- min(knots) < 0 || max(knots) > 1
  } else if (length(normalize) == 2L) {
    dim(normalize) <- c(1L, 2L)
  }
  if (is.matrix(normalize)) {
    if (M %% nrow(normalize) != 0L || ncol(normalize) != 2L) {
      stop("normalize must but be a n x 2 matrix with ", M, " divisible by n(", nrow(normalize), ")")
    }
    wa <- normalize[, 1L]
    wb <- normalize[, 2L]
    normfun <- function(x) wa * (x - wb)
    normalize <- TRUE
  } else if (normalize) {
    wa <- as.vector(1 / diff(apply(knots, 1L, range)))
    wb <- apply(knots, 1L, min)
    normfun <- function(x) wa * (x - wb)
    knots <- normfun(knots)
    normalize <- TRUE
  } else {
    normfun <- identity
    normalize <- FALSE
  }

  # trickery to get it in place
  phi <- function(x) {
    eval.parent(as.call(list(quote(.Call), C_phifunc, substitute(x), k, getOption("chebpol.threads"))))
  }

  #A <- phi(apply(knots,2, function(ck) colSums((ck-knots)^2)))
  # one day I will look into a faster solver,
  A <- phi(.Call(C_sqdiffs, knots, knots, getOption("chebpol.threads")))
  B <- rbind(1, knots)
  mat <- matrix(0, nrow(A) + nrow(B), nrow(A) + nrow(B))
  mat[seq_len(nrow(A)), seq_len(ncol(A))] <- A
  mat[(nrow(A) + 1L):nrow(mat), seq_len(ncol(B))] <- B
  mat[seq_len(ncol(B)), (ncol(A) + 1L):ncol(mat)] <- t(B)
  #    mat <- cbind(rbind(A,B),rbind(t(B),matrix(0,M+1,M+1)))
  rhs <- c(val, rep(0, M + 1L))
  wv <- try(solve(mat, rhs), silent = TRUE)
  if (inherits(wv, "try-error")) {
    if (!nowarn) {
      warning(
        "Failed to fit exactly, fallback to least squares fit.",
        if (!normalize) " You could try normalize=TRUE." else ""
      )
    }
    wv <- stats::lm.fit(mat, rhs)$coefficients
    wv[is.na(wv)] <- 0
  }
  w <- wv[seq_len(N)]
  v <- wv[(N + 1L):length(wv)]
  rm(mat, rhs, wv)
  rm(A, B)
  vectorfun(
    function(x, threads = getOption("chebpol.threads")) .Call(C_evalpolyh, normfun(x), knots, w, v, k, threads, NULL),
    arity = M,
    domain = data.frame(apply(knots, 1L, range))
  )
}
