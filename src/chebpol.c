#include "chebpol.h"
#define UNUSED(x) (void)(x)
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

// The DCT-II of x with FFTW, scaled to Chebyshev coefficients unless dct. Returns false if FFTW can't plan it,
// and always when chebpol is built without FFTW.
#ifdef HAVE_FFTW
static bool dct_fftw(const double *x, const int *dims, const int ndim, const R_xlen_t siz, double *F, bool dct,
                     int threads) {
  int rdims[ndim];
  fftw_r2r_kind kind[ndim];
  for (int i = 0; i < ndim; i++)
    kind[i] = FFTW_REDFT10; // type II DCT

  // reverse the dimensions. fftw uses row-major order.
  for (int i = 0; i < ndim; i++)
    rdims[i] = dims[ndim - i - 1];

  fftw_plan plan = fftw_plan_r2r(ndim, rdims, (double *)x, F, kind, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  if (plan != NULL) {
    fftw_execute(plan);
    fftw_destroy_plan(plan);
  } else {
    // The multi-dimensional plan can fail, e.g. with the FFTW interface of MKL. Fall back to strided 1-D transforms.
    memcpy(F, x, siz * sizeof(double));
    int stride = 1;
    for (int i = 0; i < ndim; i++) {
      fftw_iodim dim;
      fftw_plan plan;
      fftw_r2r_kind kind = FFTW_REDFT10;
      const int tlen = dims[i];
      const R_xlen_t ntrans = siz / tlen;
      dim.n = tlen;
      dim.is = stride;
      dim.os = stride;
      plan = fftw_plan_guru_r2r(1, &dim, 0, NULL, F, F, &kind, FFTW_ESTIMATE);
      if (plan == NULL) {
        if (tlen > 10000)
          Rf_warning("FFTW can't transform a dimension of length %d, falling back to a matrix product. This happens "
                     "with the FFTW interface of MKL.",
                     tlen);
        return false;
      }
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
#endif
      for (R_xlen_t j = 0; j < ntrans; j++) {
        const R_xlen_t offset = j % stride + (R_xlen_t)stride * tlen * (j / stride);
        fftw_execute_r2r(plan, F + offset, F + offset);
      }
      fftw_destroy_plan(plan);
      stride *= tlen;
    }
  }
  // adjust scale to fit our setup
  if (!dct) {
    const double isiz = 1.0 / siz;
    for (R_xlen_t i = 0; i < siz; i++)
      F[i] *= isiz;
  }
  return true;
}
#else
static bool dct_fftw(const double *x, const int *dims, const int ndim, const R_xlen_t siz, double *F, bool dct,
                     int threads) {
  return false;
}
#endif

// The same transform without FFTW, by multiplying with a cosine matrix in each dimension
static void dct_matrix(const double *x, const int *dims, const int ndim, const R_xlen_t siz, double *F, bool dct) {
  double *src, *dest, *buf;
  double **mat;
  double beta = 0;

  buf = (double *)R_alloc(siz, sizeof(double));
  mat = (double **)R_alloc(ndim, sizeof(double *));
  // Create the needed transform matrices
  // reuse same dimension matrices
  for (int j = 0; j < ndim; j++) {
    int N = dims[j];
    if (N > 10000)
      Rf_warning("A dimension of length %d needs a large matrix without FFTW, which may run out of memory.", N);
    double *jmat = NULL;
    for (int k = 0; k < j; k++) {
      if (dims[k] == N) {
        jmat = mat[k];
        break;
      }
    }
    if (jmat != NULL) {
      mat[j] = jmat;
      continue;
    }
    jmat = mat[j] = (double *)R_alloc((size_t)N * N, sizeof(double));
    for (int k = 0; k < N; k++) {
      double *jkvec = &jmat[(R_xlen_t)k * N];
      for (int i = 0; i < N; i++) {
        jkvec[i] = cospi(k * (i + 0.5) / N);
      }
    }
  }
  // We will switch src and dest between F and buf
  // We should end up with dest=F, so if
  // ndim is odd we should start with src = F, dest=buf
  if ((ndim & 1) == 1) {
    src = F;
    dest = buf;
  } else {
    src = buf;
    dest = F;
  }
  memcpy(dest, x, siz * sizeof(double));
  for (int i = ndim - 1; i >= 0; i--) {
    // transformation of dimension i, put in front
    int N = dims[i]; // length of transform
    int stride = (int)(siz / N);
    double alpha = dct ? 2.0 : 2.0 / N; // Constant for transform
    double *p = src;
    src = dest;
    dest = p;

    F77_CALL(dgemm)("TRANS", "TRANS", &N, &stride, &N, &alpha, mat[i], &N, src, &stride, &beta, dest, &N FCONE FCONE);
  }
}

static void chebyshev_coef(const double *x, const int *dims, const int ndim, double *F, bool dct, int threads) {
  R_xlen_t siz = 1;
  for (int i = 0; i < ndim; i++)
    siz *= dims[i];

  if (!dct_fftw(x, dims, ndim, siz, F, dct, threads))
    dct_matrix(x, dims, ndim, siz, F, dct);

  if (!dct) {
    // We need to adjust the first element in each dimension by 0.5
    // The DCT-II isn't exactly equal to the Chebyshev transform
    R_xlen_t blocklen = 1; // length of block
    R_xlen_t stride = 1;   // distance between blocks
    for (int i = 0; i < ndim; i++) {
      stride *= dims[i];
      for (R_xlen_t j = 0; j < siz; j += stride) {
        for (R_xlen_t s = 0; s < blocklen; s++) {
          F[j + s] *= 0.5;
        }
      }
      blocklen = stride;
    }
  }
}

SEXP c_chebyshev_coef(SEXP s_x, SEXP s_dct, SEXP s_threads) {
  int threads = Rf_asInteger(s_threads);
  bool dct = Rf_asLogical(s_dct);
  SEXP s_dim = Rf_getAttrib(s_x, R_DimSymbol);
  if (XLENGTH(s_x) > INT_MAX)
    Rf_error("`val` must have at most %d elements, not %lld.", INT_MAX, (long long)XLENGTH(s_x));
  // If no dim-attribute, assume one-dimensional
  const int len = (int)XLENGTH(s_x);
  int ndim = 1;
  const int *dims = &len;
  if (!Rf_isNull(s_dim)) {
    ndim = LENGTH(s_dim);
    if (ndim <= 0)
      Rf_error("`val` must have at least one dimension.");
    dims = INTEGER_RO(s_dim);
  }
  R_xlen_t siz = 1;
  for (int i = 0; i < ndim; i++)
    siz *= dims[i];
  if (siz == 0)
    Rf_error("`val` must have at least one element.");

  SEXP out = PROTECT(Rf_allocVector(REALSXP, siz));
  chebyshev_coef(REAL_RO(s_x), dims, ndim, REAL(out), dct, threads);
  Rf_setAttrib(out, R_DimSymbol, s_dim);
  Rf_setAttrib(out, R_DimNamesSymbol, Rf_getAttrib(s_x, R_DimNamesSymbol));
  UNPROTECT(1);

  return out;
}

static double eval_fh(const double *fv, const double *x, const double **knots, const int *dims, const int ndim,
                      const double **weights) {
  if (ndim == 0)
    return fv[0];
  R_xlen_t siz = 1;
  const int last = ndim - 1;
  const int N = dims[last];
  const double xx = x[last];
  const double *kn = knots[last];
  const double *w = weights[last];
  for (int i = 0; i < last; i++)
    siz *= dims[i];
  double num = 0, denom = 0;

  // x lies on a knot, to within the rounding of a grid of this extent
  const double tol = 10.0 * DBL_EPSILON * fabs(kn[N - 1] - kn[0]);
  for (int i = 0; i < N; i++) {
    if (fabs(xx - kn[i]) <= tol)
      return eval_fh(&fv[i * siz], x, knots, dims, last, weights);
  }

  // save a recursion step
  if (last == 0) {
    for (int i = 0; i < N; i++) {
      const double val = fv[i];
      double pole = w[i] / (xx - kn[i]);
      num += pole * val;
      denom += pole;
    }
    return num / denom;
  }
  for (int i = 0; i < N; i++) {
    const double val = eval_fh(&fv[i * siz], x, knots, dims, last, weights);
    double pole = w[i] / (xx - kn[i]);
    num += pole * val;
    denom += pole;
  }
  return num / denom;
}

SEXP c_eval_fh(SEXP s_x, SEXP s_grid, SEXP s_val, SEXP s_weights, SEXP s_threads) {
  R_xlen_t siz = 1;
  const double *val = REAL_RO(s_val);
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  SEXP s_dim = Rf_getAttrib(s_val, R_DimSymbol);
  const int *dims = INTEGER_RO(s_dim);
  const int ndim = LENGTH(s_dim);
  if (ndim <= 0)
    Rf_error("`val` must have at least one dimension.");
  if (Rf_isMatrix(s_x)) {
    if (ndim != Rf_nrows(s_x))
      Rf_error("`x` must have %d rows, one per dimension, not %d.", ndim, Rf_nrows(s_x));
  } else {
    if (ndim != LENGTH(s_x))
      Rf_error("`x` must have length %d, one per dimension, not %d.", ndim, LENGTH(s_x));
  }
  // This shouldn't happen, but we check it anyway since we'll bomb if it's wrong
  for (int i = 0; i < ndim; i++)
    siz *= dims[i];
  if (XLENGTH(s_val) != siz)
    Rf_error("`val` must have one value per grid point (%lld), not %lld.", (long long)siz, (long long)XLENGTH(s_val));

  if (LENGTH(s_grid) != ndim || LENGTH(s_weights) != ndim)
    Rf_error("`grid` and the weights must have %d dimensions, one per dimension of `val`, not %d and %d.", ndim,
             LENGTH(s_grid), LENGTH(s_weights));
  const double **knots = (const double **)R_alloc(ndim, sizeof(double *));
  const double **weights = (const double **)R_alloc(ndim, sizeof(double *));
  for (int i = 0; i < ndim; i++) {
    SEXP axis = VECTOR_ELT(s_grid, i), w = VECTOR_ELT(s_weights, i);
    if (XLENGTH(axis) != dims[i] || XLENGTH(w) != dims[i])
      Rf_error("`grid` dimension %d and its weights must have %d points, one per value along it, not %lld and %lld.",
               i + 1, dims[i], (long long)XLENGTH(axis), (long long)XLENGTH(w));
    knots[i] = REAL_RO(axis);
    weights[i] = REAL_RO(w);
  }
  const double *x = REAL_RO(s_x);
  const int n_points = Rf_isMatrix(s_x) ? Rf_ncols(s_x) : 1;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_points));
  double *p_out = REAL(out);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (n_points > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < n_points; i++) {
    p_out[i] = eval_fh(val, x + i * ndim, knots, dims, ndim, weights);
  }
  UNPROTECT(1);
  return out;
}

static void fh_weights(const int n, const double *grid, const int degree, double *w, int threads) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
#endif
  for (int k = 0; k <= n; k++) {
    const int start = (k < degree) ? 0 : k - degree, end = (k < n - degree) ? k : n - degree;
    double sum = 0.0;
    for (int i = start; i <= end; i++) {
      double prod = 1.0;
      for (int j = i; j <= i + degree; j++) {
        if (j == k)
          continue;
        prod *= grid[k] - grid[j];
      }
      sum += 1.0 / fabs(prod);
    }
    w[k] = sum * ((k % 2 != degree % 2) ? -1.0 : 1.0);
  }
}

// Compute the weights for Floater-Hormann. Formula (18) of the FH-paper
SEXP c_fh_weights(SEXP s_grid, SEXP s_degree, SEXP s_threads) {
  if (LENGTH(s_grid) != LENGTH(s_degree))
    Rf_error("`degree` must have length %d, one per dimension of `grid`, not %d.", LENGTH(s_grid), LENGTH(s_degree));

  const int ndim = LENGTH(s_grid);
  const double *grid[ndim];
  int *dims = (int *)R_alloc(ndim, sizeof(int));
  for (int i = 0; i < ndim; i++) {
    grid[i] = REAL_RO(VECTOR_ELT(s_grid, i));
    dims[i] = LENGTH(VECTOR_ELT(s_grid, i));
  }
  double *wlist[ndim];
  SEXP out = PROTECT(Rf_allocVector(VECSXP, ndim));
  for (int i = 0; i < ndim; i++) {
    SET_VECTOR_ELT(out, i, Rf_allocVector(REALSXP, dims[i]));
    wlist[i] = REAL(VECTOR_ELT(out, i));
  }

  SEXP s_degree_int = PROTECT(Rf_coerceVector(s_degree, INTSXP));
  const int *degree = INTEGER_RO(s_degree_int);
  int threads = Rf_asInteger(s_threads);
  for (int r = 0; r < ndim; r++) {
    // The paper use index from 0 to n, i.e. n+1 knots in each dimension. Therefore dims[r]-1
    fh_weights(dims[r] - 1, grid[r], degree[r], wlist[r], threads);
  }

  UNPROTECT(2);
  return out;
}

static double eval_chebyshev(const double *cf, const double *x, const int *dims, const int ndim) {
  if (ndim == 0)
    return cf[0];
  // Otherwise, use the Clenshaw algorithm
  R_xlen_t siz = 1;
  const int last = ndim - 1;
  const int N = dims[last];
  double x2 = 2.0 * x[last], bn1 = 0, bn2 = 0, bn = 0;

  // not needed for correctness, it saves a recursion level
  if (last == 0) {
    for (int i = N - 1; i > 0; i--) {
      bn2 = bn1;
      bn1 = bn;
      bn = x2 * bn1 - bn2 + cf[i];
    }
    return x[0] * bn - bn1 + cf[0];
  }
  for (int i = 0; i < last; i++)
    siz *= dims[i];
  for (int i = N - 1; i >= 0; i--) {
    bn2 = bn1;
    bn1 = bn;
    bn = x2 * bn1 - bn2 + eval_chebyshev(&cf[i * siz], x, dims, last);
  }
  return bn - x[last] * bn1;
}

SEXP c_eval_chebyshev(SEXP s_x, SEXP s_coef, SEXP s_threads) {
  if (XLENGTH(s_coef) == 0)
    Rf_error("`coef` must have at least one element.");
  R_xlen_t siz = 1;
  const double *cf = REAL_RO(s_coef);
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  SEXP s_dim = Rf_getAttrib(s_coef, R_DimSymbol);
  const int *dims = INTEGER_RO(s_dim);
  const int ndim = LENGTH(s_dim);
  if (ndim <= 0)
    Rf_error("`coef` must have at least one dimension.");
  if (Rf_isMatrix(s_x)) {
    if (ndim != Rf_nrows(s_x))
      Rf_error("`x` must have %d rows, one per dimension, not %d.", ndim, Rf_nrows(s_x));
  } else {
    if (ndim != LENGTH(s_x))
      Rf_error("`x` must have length %d, one per dimension, not %d.", ndim, LENGTH(s_x));
  }
  // This shouldn't happen, but we check it anyway since we'll bomb if it's wrong
  for (int i = 0; i < ndim; i++)
    siz *= dims[i];
  if (XLENGTH(s_coef) != siz)
    Rf_error("`coef` must have %lld elements, the product of its dimensions, not %lld.", (long long)siz,
             (long long)XLENGTH(s_coef));

  const double *x = REAL_RO(s_x);
  const int n_points = Rf_isMatrix(s_x) ? Rf_ncols(s_x) : 1;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_points));
  double *p_out = REAL(out);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (n_points > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < n_points; i++) {
    p_out[i] = eval_chebyshev(cf, x + i * ndim, dims, ndim);
  }
  UNPROTECT(1);
  return out;
}

static double eval_polyharmonic(const double *x, const double *knots, const double *weights, const double *lweights,
                                const int ndim, const int nknots, const double k) {
  // Find squared distance to each knot
  int ki = (int)k;
  double wsum = 0.0;
  if (k < 0) {
    for (R_xlen_t i = 0; i < nknots; i++) {
      const double *kn = &knots[i * ndim];
      double sqdist = 0.0;
      for (int j = 0; j < ndim; j++)
        sqdist += (x[j] - kn[j]) * (x[j] - kn[j]);
      wsum += weights[i] * exp(k * sqdist);
    }
  } else if (ki % 2 == 1) {
    for (R_xlen_t i = 0; i < nknots; i++) {
      const double *kn = &knots[i * ndim];
      double sqdist = 0.0;
      for (int j = 0; j < ndim; j++)
        sqdist += (x[j] - kn[j]) * (x[j] - kn[j]);
      if (sqdist == 0)
        continue;
      wsum += weights[i] * R_pow_di(sqrt(sqdist), ki);
    }
  } else {
    for (R_xlen_t i = 0; i < nknots; i++) {
      const double *kn = &knots[i * ndim];
      double sqdist = 0.0;
      for (int j = 0; j < ndim; j++)
        sqdist += (x[j] - kn[j]) * (x[j] - kn[j]);
      if (sqdist == 0)
        continue;
      wsum += weights[i] * log(sqrt(sqdist)) * R_pow_di(sqrt(sqdist), ki);
    }
  }
  // then the linear part, constant is the first element
  wsum += lweights[0];
  for (int i = 0; i < ndim; i++)
    wsum += lweights[i + 1] * x[i];
  return wsum;
}
SEXP c_eval_polyharmonic(SEXP s_x, SEXP s_knots, SEXP s_weights, SEXP s_lweights, SEXP s_k, SEXP s_threads) {
  const double *x = REAL_RO(s_x);
  const double *knots = REAL_RO(s_knots);
  const double *w = REAL_RO(s_weights);
  const double *lw = REAL_RO(s_lweights);
  double k = Rf_asReal(s_k);
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  const int ndim = Rf_nrows(s_knots);
  if (LENGTH(s_lweights) != ndim + 1)
    Rf_error("The linear weights must have length %d, one more than the number of dimensions, not %d.", ndim + 1,
             LENGTH(s_lweights));
  if (LENGTH(s_weights) != Rf_ncols(s_knots))
    Rf_error("The weights must have length %d, one per knot, not %d.", Rf_ncols(s_knots), LENGTH(s_weights));
  if (Rf_isMatrix(s_x)) {
    if (Rf_nrows(s_x) != Rf_nrows(s_knots))
      Rf_error("`x` must have %d rows, one per dimension, not %d.", Rf_nrows(s_knots), Rf_nrows(s_x));
  } else {
    if (Rf_length(s_x) != Rf_nrows(s_knots))
      Rf_error("`x` must have length %d, one per dimension, not %d.", Rf_nrows(s_knots), Rf_length(s_x));
  }
  const int n_points = Rf_isMatrix(s_x) ? Rf_ncols(s_x) : 1;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_points));
  double *p_out = REAL(out);
  const int nknots = Rf_ncols(s_knots);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (n_points > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < n_points; i++) {
    p_out[i] = eval_polyharmonic(x + i * ndim, knots, w, lw, ndim, nknots, k);
  }
  UNPROTECT(1);
  return out;
}

// an R-free evalongrid. Recursive
static void eval_on_grid(void (*fun)(double *x, double *y, int valuedim, void *ud), double *arg, const double **grid,
                         const int *dims, const int ndim, const int valuedim, double *result, void *userdata) {
  int last = ndim - 1;
  R_xlen_t stride = valuedim;

  if (ndim == 0) {
    fun(arg, &result[0], valuedim, userdata);
    return;
  }

  for (int i = 0; i < last; i++)
    stride *= dims[i];
  for (int i = 0; i < dims[last]; i++) {
    arg[last] = grid[last][i];
    eval_on_grid(fun, arg, grid, dims, last, valuedim, &result[i * stride], userdata);
  }
}

// The call to the R function, and its value at the first grid point until the traversal has used it
typedef struct {
  SEXP call;
  SEXP first;
} grid_call;

static void call_fun(double *x, double *y, const int valuedim, void *userdata) {
  grid_call *gc = (grid_call *)userdata;
  // The first grid point was evaluated to find the length of the value, don't call the function there again
  if (gc->first != NULL) {
    memcpy(y, REAL_RO(gc->first), valuedim * sizeof(double));
    gc->first = NULL;
    return;
  }
  SEXP call = gc->call;
  // Reuse the argument vector of the previous call, unless the function kept a reference to it.
  // The call holds one reference itself, so the vector is shared only if something else does too.
  SEXP arg = CADR(call);
  const R_xlen_t ndim = XLENGTH(arg);
  if (MAYBE_SHARED(arg)) {
    arg = Rf_allocVector(REALSXP, ndim);
    SETCADR(call, arg);
  }
  if (ndim > 0)
    memcpy(REAL(arg), x, ndim * sizeof(double));
  PROTECT_INDEX ipx;
  SEXP fv = Rf_eval(call, R_BaseEnv);
  PROTECT_WITH_INDEX(fv, &ipx);
  if (!Rf_isNumeric(fv))
    Rf_error("`fun` must return a numeric value.");
  if (LENGTH(fv) != valuedim)
    Rf_error("`fun` must return a value of length %d at every grid point, not %d.", valuedim, LENGTH(fv));
  REPROTECT(fv = Rf_coerceVector(fv, REALSXP), ipx);
  memcpy(y, REAL_RO(fv), valuedim * sizeof(double));
  UNPROTECT(1);
}

SEXP c_eval_on_grid(SEXP s_fun, SEXP s_grid) {
  if (!Rf_isFunction(s_fun))
    Rf_error("`fun` must be a function.");
  const int ndim = LENGTH(s_grid);
  if (ndim == 0)
    Rf_error("`grid` must have at least one dimension.");
  const double *grid[ndim];
  R_xlen_t len = 1;
  SEXP s_dim = PROTECT(Rf_allocVector(INTSXP, ndim));
  int *dims = INTEGER(s_dim);
  for (int i = 0; i < ndim; i++) {
    grid[i] = REAL_RO(VECTOR_ELT(s_grid, i));
    dims[i] = LENGTH(VECTOR_ELT(s_grid, i));
    if (dims[i] == 0)
      Rf_error("`grid` dimension %d must have at least one point.", i + 1);
    if (len > R_XLEN_T_MAX / dims[i])
      Rf_error("`grid` has too many points.");
    len *= dims[i];
  }

  // Create a call to s_fun
  SEXP s_call = PROTECT(Rf_lang2(s_fun, R_NilValue));
  // and an argument list to it
  SEXP s_arg = PROTECT(Rf_allocVector(REALSXP, ndim));
  SETCADR(s_call, s_arg);

  // find the dimension of the result from the value at the first grid point
  for (int i = 0; i < ndim; i++)
    REAL(s_arg)[i] = grid[i][0];
  SEXP fv = PROTECT(Rf_eval(s_call, R_BaseEnv));
  if (!Rf_isNumeric(fv))
    Rf_error("`fun` must return a numeric value.");
  const int valuedim = LENGTH(fv);
  grid_call gc = {s_call, PROTECT(Rf_coerceVector(fv, REALSXP))};
  if (valuedim > 0 && len > R_XLEN_T_MAX / valuedim)
    Rf_error("`fun` returns too many values on `grid` to store them.");

  // Now, fill in the grid points one by one in a buffer, call_fun copies it into the call
  SEXP out = PROTECT(Rf_allocVector(REALSXP, len * valuedim));
  double *arg = (double *)R_alloc(ndim, sizeof(double));
  eval_on_grid(call_fun, arg, grid, dims, ndim, valuedim, REAL(out), (void *)&gc);

  // Set the dim-attribute
  if (valuedim == 1) {
    Rf_setAttrib(out, R_DimSymbol, s_dim);
  } else {
    SEXP s_new_dim = PROTECT(Rf_allocVector(INTSXP, ndim + 1));
    int *newdim = INTEGER(s_new_dim);
    newdim[0] = valuedim;
    memcpy(newdim + 1, dims, ndim * sizeof(int));
    Rf_setAttrib(out, R_DimSymbol, s_new_dim);
    UNPROTECT(1);
  }
  UNPROTECT(6);
  return out;
}

static double eval_multilinear(const int ndim, const double *x, const double **grid, const int *dims,
                               const double *values, const blend_type blend) {

  // the blended weights of the lower and upper end of the interval in each dimension
  double wlo[ndim], wup[ndim];
  R_xlen_t valpos = 0;
  double ipval = 0;
  // For each dimension, find the grid interval holding x and the weight within it, and accumulate the
  // linear position valpos of the upper corner. Points outside the grid use the outermost interval.

  R_xlen_t stride = 1;
  for (int i = 0; i < ndim; i++) {
    if (ISNAN(x[i]))
      return x[i];
    const double *gdvec = grid[i];
    // upper end of the grid interval holding x[i]
    const int gp = find_interval(gdvec, dims[i], x[i]) + 1;
    const double weight = (x[i] - gdvec[gp - 1]) / (gdvec[gp] - gdvec[gp - 1]);
    wup[i] = blend_fun(weight, blend);
    wlo[i] = blend_fun(1 - weight, blend);
    valpos += stride * gp;
    stride *= dims[i];
  }

  // loop over the corners of the box, sum values with weights
  double wsum = 0.0;
  for (int i = 0; i < (1 << ndim); i++) {
    // i represents a corner. bit=1 if upper corner, 0 if lower corner.
    // We should find its weight
    // it's a product of the weights from each dimension
    R_xlen_t vpos = valpos;
    R_xlen_t stride = 1;
    double cw = 1;
    for (int g = 0; g < ndim; g++) {
      if ((1 << g) & i) {
        cw *= wup[g];
      } else {
        cw *= wlo[g];
        vpos -= stride;
      }
      stride *= dims[g];
    }
    wsum += cw;
    ipval += cw * values[vpos];
  }
  return ipval / wsum;
}

SEXP c_eval_multilinear(SEXP s_x, SEXP s_grid, SEXP s_val, SEXP s_blend, SEXP s_threads) {
  const int ndim = LENGTH(s_grid);
  R_xlen_t gridsize = 1;
  int dims[ndim];
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  const double *grid[ndim];
  blend_type blend = BLEND_LINEAR;
  if (!Rf_isReal(s_val))
    Rf_error("`val` must be numeric.");
  if (!Rf_isReal(s_x))
    Rf_error("`x` must be numeric.");
  if (Rf_isMatrix(s_x) && Rf_nrows(s_x) != ndim)
    Rf_error("`x` must have %d rows, one per dimension, not %d.", ndim, Rf_nrows(s_x));
  if (!Rf_isMatrix(s_x) && LENGTH(s_x) != ndim)
    Rf_error("`x` must have length %d, one per dimension, not %d.", ndim, LENGTH(s_x));

  for (int i = 0; i < ndim; i++) {
    dims[i] = LENGTH(VECTOR_ELT(s_grid, i));
    grid[i] = REAL_RO(VECTOR_ELT(s_grid, i));
    gridsize *= dims[i];
  }

  if (XLENGTH(s_val) != gridsize)
    Rf_error("`val` must have one value per grid point (%lld), not %lld.", (long long)gridsize,
             (long long)XLENGTH(s_val));
  if (!Rf_isNull(s_blend))
    blend = (blend_type)Rf_asInteger(s_blend);
  const int n_points = Rf_isMatrix(s_x) ? Rf_ncols(s_x) : 1;
  const double *x = REAL_RO(s_x);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_points));
  double *p_out = REAL(out);
  const double *val = REAL_RO(s_val);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (n_points > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < n_points; i++) {
    p_out[i] = eval_multilinear(ndim, x + i * ndim, grid, dims, val, blend);
  }
  UNPROTECT(1);
  return out;
}

SEXP c_sq_dists(SEXP s_x1, SEXP s_x2, SEXP s_threads) {
  // The squared distance from each column of s_x1 to each column of s_x2, with one row per column of s_x1
  const int r1 = Rf_nrows(s_x1), c1 = Rf_ncols(s_x1), r2 = Rf_nrows(s_x2), c2 = Rf_ncols(s_x2);
  const R_xlen_t N = (R_xlen_t)c1 * c2;
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  SEXP out = PROTECT(Rf_allocVector(REALSXP, N));
  double *p_out = REAL(out);
  const double *x1 = REAL_RO(s_x1), *x2 = REAL_RO(s_x2);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1 && c2 > 1)
#endif
  for (R_xlen_t j = 0; j < c2; j++) {
    const double *x2p = x2 + j * r2;
    for (R_xlen_t i = 0; i < c1; i++) {
      const double *x1p = x1 + i * r1;
      double csum = 0.0;
      for (int k = 0; k < r1; k++) {
        csum += (x1p[k] - x2p[k]) * (x1p[k] - x2p[k]);
      }
      p_out[i + j * c1] = csum;
    }
  }
  SEXP s_new_dim = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(s_new_dim)[0] = c1;
  INTEGER(s_new_dim)[1] = c2;
  Rf_setAttrib(out, R_DimSymbol, s_new_dim);
  UNPROTECT(2);
  return out;
}

// A uniform grid over the knots. Cell i lists the simplices items[offsets[i]] .. items[offsets[i + 1] - 1]
// whose bounding box overlaps it.
typedef struct {
  const double *lo, *hi, *inv;
  int m;
  const int *offsets, *items;
} simplex_grid;

static int grid_cell(const double v, const double lo, const double inv, const int m) {
  const int c = (int)((v - lo) * inv);
  return c < 0 ? 0 : (c >= m ? m - 1 : c);
}

// Whether simplex holds x. If so, vec has the barycentric coordinates of x.
static bool in_simplex(const double *x, const R_xlen_t simplex, const double *lumats, const int *ipivs,
                       const double *bbox, const int ndim, double *vec) {
  // x must be in the bounding box
  const double *box = bbox + simplex * 2 * ndim;
  for (int i = 0; i < ndim; i++) {
    if (x[i] < box[0] || x[i] > box[1])
      return false;
    box += 2;
  }
  // bounding box matches. Transform to barycentric coordinates, check that they are positive.
  int N = ndim + 1, one = 1, info;
  memcpy(vec, x, ndim * sizeof(double));
  vec[ndim] = 1.0;
  F77_CALL(dgetrs)("N", &N, &one, lumats + simplex * N * N, &N, ipivs + simplex * N, vec, &N, &info FCONE);
  for (int d = 0; d < N; d++)
    if (vec[d] < -1e-10)
      return false;
  return true;
}

// The blended value in a simplex with vertices tri, at barycentric coordinates vec
static double simplex_value(const double *vec, const int *tri, const double *val, const int N, const blend_type blend) {
  if (blend == BLEND_SQUARE) {
    // the value at the vertex with the largest barycentric coordinate
    int largest = 0;
    for (int d = 1; d < N; d++)
      if (vec[d] > vec[largest])
        largest = d;
    return val[tri[largest] - 1];
  }
  double sum = 0;
  double wsum = 0;
  for (int d = 0; d < N; d++) {
    // rounding can put a barycentric coordinate slightly outside [0, 1]; <= also turns -0.0 into 0
    const double coord = vec[d] <= 0 ? 0 : (vec[d] >= 1 ? 1 : vec[d]);
    const double w = blend_fun(coord, blend);
    if (w == 0)
      continue;
    wsum += w;
    sum += w * val[tri[d] - 1];
  }
  return sum / wsum;
}

static double eval_simplexlinear(const double *x, const double *knots, const int *dtri, const double *lumats,
                                 const int *ipivs, const double *bbox, bool epol, const double *val, const int ndim,
                                 const int numsimplex, int nknots, blend_type blend, const simplex_grid *sg) {

  for (int i = 0; i < ndim; i++)
    if (ISNAN(x[i]))
      return x[i];

  double vec[ndim + 1];
  int N = ndim + 1, one = 1, info;
  // Only the simplices whose bounding box overlaps the grid cell of x can hold it. They are listed in
  // increasing order, so the first hit is the same simplex a scan of all of them finds.
  const int *cand = NULL;
  int ncand = numsimplex;
  if (sg != NULL) {
    R_xlen_t cell = 0, stride = 1;
    for (int i = 0; i < ndim; i++) {
      if (x[i] < sg->lo[i] || x[i] > sg->hi[i]) {
        cell = -1;
        break;
      }
      cell += stride * grid_cell(x[i], sg->lo[i], sg->inv[i], sg->m);
      stride *= sg->m;
    }
    cand = cell < 0 ? NULL : sg->items + sg->offsets[cell];
    ncand = cell < 0 ? 0 : sg->offsets[cell + 1] - sg->offsets[cell];
  }
  for (int c = 0; c < ncand; c++) {
    const R_xlen_t simplex = cand != NULL ? cand[c] : c;
    if (in_simplex(x, simplex, lumats, ipivs, bbox, ndim, vec))
      return simplex_value(vec, dtri + simplex * N, val, N, blend);
  }

  if (epol) {
    // Find the closest knot
    double mindist = INFINITY;
    int nearknot = 0;
    for (int k = 0; k < nknots; k++) {
      double dist = 0.0;
      for (int i = 0; i < ndim; i++)
        dist += (x[i] - knots[(R_xlen_t)k * ndim + i]) * (x[i] - knots[(R_xlen_t)k * ndim + i]);
      if (dist < mindist) {
        mindist = dist;
        nearknot = k;
      }
    }

    R_xlen_t nearest = -1;
    for (R_xlen_t simplex = 0; simplex < numsimplex; simplex++) {
      for (int j = 0; j < N; j++)
        if (dtri[simplex * N + j] == nearknot + 1) {
          nearest = simplex;
          break;
        }
      if (nearest != -1)
        break;
    }

    if (nearest == -1)
      return NA_REAL;
    const double *lumat = lumats + nearest * N * N;
    const int *ipiv = ipivs + nearest * N;
    const int *tri = dtri + nearest * (ndim + 1);
    memcpy(vec, x, ndim * sizeof(double));
    vec[ndim] = 1.0;
    F77_CALL(dgetrs)("N", &N, &one, lumat, &N, ipiv, vec, &N, &info FCONE);
    double sum = 0;
    for (int d = 0; d <= ndim; d++)
      sum += vec[d] * val[tri[d] - 1];
    return sum;
  }
  return NA_REAL;
}

SEXP c_eval_simplexlinear(SEXP s_x, SEXP s_knots, SEXP s_dtri, SEXP s_adata, SEXP s_val, SEXP s_extrapolate,
                          SEXP s_blend, SEXP s_threads) {
  const int *dtri = INTEGER_RO(s_dtri);
  const int numsimplex = Rf_ncols(s_dtri);
  const double *lumats = REAL_RO(VECTOR_ELT(s_adata, 0));
  const int *pivots = INTEGER_RO(VECTOR_ELT(s_adata, 1));
  const double *bboxmat = REAL_RO(VECTOR_ELT(s_adata, 2));
  const double *x = REAL_RO(s_x);
  const int ndim = Rf_nrows(s_knots);
  const int n_points = Rf_ncols(s_x);
  if (Rf_nrows(s_x) != ndim)
    Rf_error("`x` must have %d rows, one per dimension, not %d.", ndim, Rf_nrows(s_x));
  const int nknots = Rf_ncols(s_knots);
  if (XLENGTH(s_val) != nknots)
    Rf_error("`val` must have one value per knot (%d), not %lld.", nknots, (long long)XLENGTH(s_val));
  const double *knots = REAL_RO(s_knots);
  const double *val = REAL_RO(s_val);
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  SEXP out = PROTECT(Rf_allocVector(REALSXP, Rf_ncols(s_x)));
  double *p_out = REAL(out);
  bool epol = Rf_asLogical(s_extrapolate);
  blend_type blend = BLEND_LINEAR;
  if (!Rf_isNull(s_blend))
    blend = (blend_type)Rf_asInteger(s_blend);
  // interpolants saved by earlier versions have no grid and scan all simplices
  simplex_grid grid, *sg = NULL;
  if (LENGTH(s_adata) >= 7) {
    const double *bounds = REAL_RO(VECTOR_ELT(s_adata, 3));
    grid.lo = bounds;
    grid.hi = bounds + ndim;
    grid.inv = bounds + 2 * ndim;
    grid.m = Rf_asInteger(VECTOR_ELT(s_adata, 4));
    grid.offsets = INTEGER_RO(VECTOR_ELT(s_adata, 5));
    grid.items = INTEGER_RO(VECTOR_ELT(s_adata, 6));
    sg = &grid;
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(guided) if (threads > 1 && n_points > 1)
#endif
  for (R_xlen_t i = 0; i < n_points; i++) {
    p_out[i] = eval_simplexlinear(x + i * ndim, knots, dtri, lumats, pivots, bboxmat, epol, val, ndim, numsimplex,
                                  nknots, blend, sg);
  }
  UNPROTECT(1);
  return out;
}

// The cells c0[j] .. c1[j] of the grid that a bounding box overlaps in each dimension, and how many cells that is
static R_xlen_t box_cells(const double *box, const double *lo, const double *inv, const int m, const int ndim, int *c0,
                          int *c1) {
  R_xlen_t ncells = 1;
  for (int j = 0; j < ndim; j++) {
    c0[j] = grid_cell(box[2 * j], lo[j], inv[j], m);
    c1[j] = grid_cell(box[2 * j + 1], lo[j], inv[j], m);
    ncells *= c1[j] - c0[j] + 1;
  }
  return ncells;
}

// The position of the cell with coordinates cc in the list of cells
static R_xlen_t cell_index(const int *cc, const int m, const int ndim) {
  R_xlen_t cell = 0, stride = 1;
  for (int j = 0; j < ndim; j++) {
    cell += stride * cc[j];
    stride *= m;
  }
  return cell;
}

// Step cc to the next cell in c0 .. c1. Returns false after the last one.
static bool next_cell(int *cc, const int *c0, const int *c1, const int ndim) {
  for (int j = 0; j < ndim; j++) {
    if (cc[j] < c1[j]) {
      cc[j]++;
      return true;
    }
    cc[j] = c0[j];
  }
  return false;
}

// Index the simplices by a uniform grid over the knots. Fills in elements 3 to 6 of out: the bounds, the
// cells per dimension, and the offsets and items of the cell lists.
static void build_simplex_grid(SEXP out, const double *bboxmat, const int ndim, const int numsimplex) {
  // the bounds of the grid, and the number of cells per unit length in each dimension
  SEXP s_bounds = Rf_allocVector(REALSXP, 3 * ndim);
  SET_VECTOR_ELT(out, 3, s_bounds);
  double *lo = REAL(s_bounds), *hi = lo + ndim, *inv = lo + 2 * ndim;
  for (int j = 0; j < ndim; j++) {
    lo[j] = INFINITY;
    hi[j] = -INFINITY;
  }
  for (R_xlen_t s = 0; s < numsimplex; s++) {
    const double *box = bboxmat + s * 2 * ndim;
    for (int j = 0; j < ndim; j++) {
      if (box[2 * j] < lo[j])
        lo[j] = box[2 * j];
      if (box[2 * j + 1] > hi[j])
        hi[j] = box[2 * j + 1];
    }
  }

  // About one cell per simplex, with m cells in each dimension, halved while the simplices overlap too many
  int *c0 = (int *)R_alloc(ndim, sizeof(int)), *c1 = (int *)R_alloc(ndim, sizeof(int));
  int *cc = (int *)R_alloc(ndim, sizeof(int));
  // the lists of all cells together hold at most 16 entries per simplex
  const R_xlen_t maxitems = numsimplex < INT_MAX / 16 ? 16 * (R_xlen_t)numsimplex : INT_MAX;
  int m = numsimplex > 0 ? (int)floor(pow((double)numsimplex, 1.0 / ndim)) : 1;
  if (m < 1)
    m = 1;
  R_xlen_t nitems;
  for (;; m /= 2) {
    for (int j = 0; j < ndim; j++)
      inv[j] = hi[j] > lo[j] ? m / (hi[j] - lo[j]) : 0.0;
    nitems = 0;
    for (R_xlen_t s = 0; s < numsimplex && nitems <= maxitems; s++)
      nitems += box_cells(bboxmat + s * 2 * ndim, lo, inv, m, ndim, c0, c1);
    if (m == 1 || nitems <= maxitems)
      break;
  }
  R_xlen_t ncells = 1;
  for (int j = 0; j < ndim; j++)
    ncells *= m;
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(m));

  SEXP s_offsets = Rf_allocVector(INTSXP, ncells + 1);
  SET_VECTOR_ELT(out, 5, s_offsets);
  SEXP s_items = Rf_allocVector(INTSXP, nitems);
  SET_VECTOR_ELT(out, 6, s_items);
  int *offsets = INTEGER(s_offsets), *items = INTEGER(s_items);

  // count the simplices of each cell
  memset(offsets, 0, (ncells + 1) * sizeof(int));
  for (R_xlen_t s = 0; s < numsimplex; s++) {
    box_cells(bboxmat + s * 2 * ndim, lo, inv, m, ndim, c0, c1);
    memcpy(cc, c0, ndim * sizeof(int));
    do {
      offsets[cell_index(cc, m, ndim) + 1]++;
    } while (next_cell(cc, c0, c1, ndim));
  }
  for (R_xlen_t i = 0; i < ncells; i++)
    offsets[i + 1] += offsets[i];

  // fill in the simplices in increasing order, next[i] is the free position in the list of cell i
  int *next = (int *)R_alloc(ncells, sizeof(int));
  memcpy(next, offsets, ncells * sizeof(int));
  for (R_xlen_t s = 0; s < numsimplex; s++) {
    box_cells(bboxmat + s * 2 * ndim, lo, inv, m, ndim, c0, c1);
    memcpy(cc, c0, ndim * sizeof(int));
    do {
      items[next[cell_index(cc, m, ndim)]++] = (int)s;
    } while (next_cell(cc, c0, c1, ndim));
  }
}

SEXP c_analyze_simplex(SEXP s_dtri, SEXP s_knots, SEXP s_threads) {
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  const int *dtri = INTEGER_RO(s_dtri);
  const double *knots = REAL_RO(s_knots);
  const int ndim = Rf_nrows(s_knots);
  const int numsimplex = Rf_ncols(s_dtri);
  int N = ndim + 1;
  // The LU factors, pivots and bounding boxes of the simplices, then the grid that indexes them
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  // Allocate everything before parallel for
  SET_VECTOR_ELT(out, 0, Rf_allocVector(REALSXP, (R_xlen_t)numsimplex * N * N));
  SET_VECTOR_ELT(out, 1, Rf_allocVector(INTSXP, (R_xlen_t)numsimplex * N));
  SET_VECTOR_ELT(out, 2, Rf_allocMatrix(REALSXP, 2 * ndim, numsimplex));
  double *lumats = REAL(VECTOR_ELT(out, 0));
  int *pivots = INTEGER(VECTOR_ELT(out, 1));
  double *bboxmat = REAL(VECTOR_ELT(out, 2));

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1 && numsimplex > 1)
#endif
  for (R_xlen_t simplex = 0; simplex < numsimplex; simplex++) {
    const int *tri = dtri + simplex * (ndim + 1);
    int *ipiv = pivots + simplex * N;
    for (int j = 0; j < N; j++)
      ipiv[j] = j + 1;
    double *lumat = lumats + simplex * N * N;
    int info;
    double *m = bboxmat + simplex * 2 * ndim;
    for (int j = 0; j < ndim; j++) {
      m[2 * j] = INFINITY;
      m[2 * j + 1] = -INFINITY;
    }
    // Copy all vertices into mat, with a final row of ones to LU-factorize
    // Preparation for transforming into barycentric coordinates
    for (int v = 0; v <= ndim; v++) {
      const double *knot = knots + (tri[v] - 1) * ndim;
      for (int j = 0; j < ndim; j++) {
        lumat[v * (ndim + 1) + j] = knot[j];
        if (knot[j] < m[2 * j])
          m[2 * j] = knot[j];
        if (knot[j] > m[2 * j + 1])
          m[2 * j + 1] = knot[j];
      }
      lumat[v * (ndim + 1) + ndim] = 1.0;
    }
    F77_CALL(dgetrf)(&N, &N, lumat, &N, ipiv, &info);
  }

  build_simplex_grid(out, bboxmat, ndim, numsimplex);

  UNPROTECT(1);
  return out;
}

SEXP c_has_fftw(void) {
#ifdef HAVE_FFTW
  return Rf_ScalarLogical(TRUE);
#else
  return Rf_ScalarLogical(FALSE);
#endif
}

// inplace to save memory, do repeated multiplication instead of pow(). It's faster.
SEXP c_rbf_kernel(SEXP s_x, SEXP s_k, SEXP s_threads) {
  double k = Rf_asReal(s_k);
  double *x = REAL(s_x);
  const R_xlen_t n = XLENGTH(s_x);
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  double *y;
  SEXP out;
  if (MAYBE_REFERENCED(s_x)) {
    out = PROTECT(Rf_allocVector(REALSXP, n));
    Rf_setAttrib(out, R_DimSymbol, Rf_getAttrib(s_x, R_DimSymbol));
    y = REAL(out);
  } else {
    out = s_x;
    y = x;
  }

  if (k < 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
#endif
    for (R_xlen_t i = 0; i < n; i++)
      y[i] = exp(k * x[i]);
  } else {
    int ki = Rf_asInteger(s_k);
    if (ki % 2 == 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
#endif
      for (R_xlen_t i = 0; i < n; i++) {
        // it's the sqrt(x) to ki'th power
        if (x[i] <= 0.0) {
          y[i] = 0.0;
          continue;
        }
        y[i] = R_pow_di(sqrt(x[i]), ki);
      }
    } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
#endif
      for (R_xlen_t i = 0; i < n; i++) {
        // it's sqrt(x) to ki'th power, multiplied by 0.5 log(x)
        if (x[i] <= 0.0) {
          y[i] = 0.0;
          continue;
        }
        y[i] = 0.5 * log(x[i]) * R_pow_di(sqrt(x[i]), ki);
      }
    }
  }
  if (y != x)
    UNPROTECT(1);
  return out;
}
#ifndef HAVE_ALGLIB
SEXP c_make_crbf(SEXP s_data, SEXP s_layers, SEXP s_rbase, SEXP s_lambda) {
  UNUSED(s_data);
  UNUSED(s_layers);
  UNUSED(s_rbase);
  UNUSED(s_lambda);
  Rf_error("method = \"crbf\" is not available because chebpol has no ALGLIB support.");
  return R_NilValue;
}
SEXP c_eval_crbf(SEXP s_x, SEXP s_model, SEXP s_threads) {
  UNUSED(s_x);
  UNUSED(s_model);
  UNUSED(s_threads);
  Rf_error("method = \"crbf\" is not available because chebpol has no ALGLIB support.");
  return R_NilValue;
}
SEXP c_has_alglib(void) { return Rf_ScalarLogical(FALSE); }

#endif
