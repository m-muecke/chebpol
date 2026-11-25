#include "chebpol.h"
#define UNUSED(x) (void)(x)
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

static void chebcoef(const double *x, const int *dims, const int rank, double *F, bool dct, int threads) {
  R_xlen_t siz = 1;
  for (int i = 0; i < rank; i++)
    siz *= dims[i];

#ifdef HAVE_FFTW
  bool brokenfftw = false;
  double isiz = 1.0 / siz;
  int rdims[rank];
  // Create a plan
  fftw_r2r_kind kind[rank];
  for (int i = 0; i < rank; i++)
    kind[i] = FFTW_REDFT10; // type II DCT

  // reverse the dimensions. fftw uses row-major order.
  for (int i = 0; i < rank; i++)
    rdims[i] = dims[rank - i - 1];

  // Plan and execute

  fftw_plan plan = fftw_plan_r2r(rank, rdims, (double *)x, F, kind, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  if (plan != NULL) {
    fftw_execute(plan);
    fftw_destroy_plan(plan);
  } else {
    // That's funny, perhaps we've been given MKL with its amputated FFTW?
    // Let's try a series of strided 1d-transforms
    memcpy(F, x, siz * sizeof(double));
    int stride = 1;
    for (int i = 0; i < rank; i++) {
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
          Rf_warning("FFTW fails on strided long (%d) vector, is it MKL?", tlen);
        brokenfftw = true;
        break;
      }
      // Do them
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
  if (!dct)
    for (R_xlen_t i = 0; i < siz; i++)
      F[i] *= isiz;

  if (brokenfftw) {
#endif
    double *src, *dest, *buf;
    double **mat;
    double beta = 0;

    // Some work space
    //
    buf = (double *)R_alloc(siz, sizeof(double));
    mat = (double **)R_alloc(rank, sizeof(double *));
    // Create the needed transform matrices
    // reuse same dimension matrices
    for (int j = 0; j < rank; j++) {
      // Is it there already?
      int N = dims[j];
      if (N > 10000)
        Rf_warning("Long vector (%d), no FFTW. This may go wrong.", N);
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
    // rank is odd we should start with src = F, dest=buf
    if ((rank & 1) == 1) {
      src = F;
      dest = buf;
    } else {
      src = buf;
      dest = F;
    }
    memcpy(dest, x, siz * sizeof(double));
    for (int i = rank - 1; i >= 0; i--) {
      // transformation of dimension i, put in front
      int N = dims[i]; // length of transform
      int stride = siz / N;
      double alpha = dct ? 2.0 : 2.0 / N; // Constant for transform
      // swap src and dest
      double *p = src;
      src = dest;
      dest = p;

      F77_CALL(dgemm)("TRANS", "TRANS", &N, &stride, &N, &alpha, mat[i], &N, src, &stride, &beta, dest, &N FCONE FCONE);
    }
#ifdef HAVE_FFTW
  }
#endif
  if (!dct) {
    // We need to adjust the first element in each dimension by 0.5
    // The DCT-II isn't exactly equal to the Chebyshev transform
    R_xlen_t blocklen = 1; // length of block
    R_xlen_t stride = 1;   // distance between blocks
    for (int i = 0; i < rank; i++) {
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

SEXP R_chebcoef(SEXP x, SEXP sdct, SEXP Sthreads) {
  int threads = Rf_asInteger(Sthreads);
  if (!Rf_isLogical(sdct) || LENGTH(sdct) < 1)
    Rf_error("dct must be a logical");
  bool dct = LOGICAL_RO(sdct)[0];
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  // If no dim-attribute, assume one-dimensional
  const int sdims = LENGTH(x);
  int rank = 1;
  const int *dims = &sdims;
  if (!Rf_isNull(dim)) {
    rank = LENGTH(dim);
    if (rank <= 0)
      Rf_error("Rank must be positive");
    dims = INTEGER_RO(dim);
  }
  R_xlen_t siz = 1;
  for (int i = 0; i < rank; i++)
    siz *= dims[i];
  if (siz == 0)
    Rf_error("array size must be positive");

  SEXP resvec = PROTECT(Rf_allocVector(REALSXP, siz));
  chebcoef(REAL_RO(x), dims, rank, REAL(resvec), dct, threads);
  Rf_setAttrib(resvec, R_DimSymbol, dim);
  Rf_setAttrib(resvec, R_DimNamesSymbol, Rf_getAttrib(x, R_DimNamesSymbol));
  UNPROTECT(1);

  return resvec;
}

static double FH(const double *fv, const double *x, const double **knots, const int *dims, const int rank,
                 const double **weights) {
  // Use Floater-Hormann method
  if (rank == 0)
    return fv[0];
  R_xlen_t siz = 1;
  const int newrank = rank - 1;
  const int N = dims[newrank];
  const double xx = x[newrank];
  const double *kn = knots[newrank];
  const double *w = weights[newrank];
  for (int i = 0; i < newrank; i++)
    siz *= dims[i];
  double num = 0, denom = 0;

  // Special case:
  for (int i = 0; i < N; i++) {
    if (fabs(xx - kn[i]) < 10.0 * DBL_EPSILON)
      return FH(&fv[i * siz], x, knots, dims, newrank, weights);
  }

  // save a recursion step
  if (newrank == 0) {
    for (int i = 0; i < N; i++) {
      const double val = fv[i];
      double pole = w[i] / (xx - kn[i]);
      num += pole * val;
      denom += pole;
    }
    return num / denom;
  }
  for (int i = 0; i < N; i++) {
    const double val = FH(&fv[i * siz], x, knots, dims, newrank, weights);
    double pole = w[i] / (xx - kn[i]);
    num += pole * val;
    denom += pole;
  }
  return num / denom;
}

SEXP R_FH(SEXP inx, SEXP vals, SEXP grid, SEXP Sweights, SEXP Rthreads, SEXP spare) {
  R_xlen_t siz = 1;
  const double *val = REAL_RO(vals);
#ifdef _OPENMP
  int threads = Rf_asInteger(Rthreads);
#endif
  UNUSED(spare);
  SEXP dim = Rf_getAttrib(vals, R_DimSymbol);
  const int *dims = INTEGER_RO(dim);
  const int rank = LENGTH(dim);
  if (rank <= 0)
    Rf_error("rank must be positive");
  if (Rf_isMatrix(inx)) {
    if (rank != Rf_nrows(inx))
      Rf_error("coefficient rank(%d) must match number of rows(%d)", LENGTH(dim), Rf_nrows(inx));
  } else {
    if (rank != LENGTH(inx))
      Rf_error("coefficient rank(%d) does not match argument length(%d)", LENGTH(dim), LENGTH(inx));
  }
  // This shouldn't happen, but we check it anyway since we'll bomb if it's wrong
  for (int i = 0; i < rank; i++)
    siz *= dims[i];
  if (XLENGTH(vals) != siz)
    Rf_error("coefficient length(%lld) does not match data length(%lld)", (long long)XLENGTH(vals), (long long)siz);

  if (LENGTH(grid) != rank)
    Rf_error("There must be one value for each knot");
  const double **knots = (const double **)R_alloc(rank, sizeof(double *));
  const double **weights = (const double **)R_alloc(rank, sizeof(double *));
  for (int i = 0; i < rank; i++) {
    knots[i] = REAL_RO(VECTOR_ELT(grid, i));
    weights[i] = REAL_RO(VECTOR_ELT(Sweights, i));
  }
  const double *xp = REAL_RO(inx);
  const int numvec = Rf_isMatrix(inx) ? Rf_ncols(inx) : 1;
  SEXP resvec = PROTECT(Rf_allocVector(REALSXP, numvec));
  double *out = REAL(resvec);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < numvec; i++) {
    out[i] = FH(val, xp + i * rank, knots, dims, rank, weights);
  }
  UNPROTECT(1);
  return resvec;
}

static inline void fhweights(const int n, const double *grid, const int d, double *w, int threads) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
#endif
  for (int k = 0; k <= n; k++) {
    const int start = (k < d) ? 0 : k - d, end = (k < n - d) ? k : n - d;
    double sum = 0.0;
    for (int i = start; i <= end; i++) {
      double prod = 1.0;
      for (int j = i; j <= i + d; j++) {
        if (j == k)
          continue;
        prod *= grid[k] - grid[j];
      }
      sum += 1.0 / fabs(prod);
    }
    w[k] = sum * ((k % 2 != d % 2) ? -1.0 : 1.0);
  }
}

// Compute the weights for Floater-Hormann. Formula (18) of the FH-paper
SEXP R_fhweights(SEXP Sgrid, SEXP Sd, SEXP Sthreads) {
  if (LENGTH(Sgrid) != LENGTH(Sd))
    Rf_error("Length of grid (%d) should equal length of d (%d)", LENGTH(Sgrid), LENGTH(Sd));

  const int rank = LENGTH(Sgrid);
  const double *grid[rank];
  int *dims = (int *)R_alloc(rank, sizeof(int));
  for (int i = 0; i < rank; i++) {
    grid[i] = REAL_RO(VECTOR_ELT(Sgrid, i));
    dims[i] = LENGTH(VECTOR_ELT(Sgrid, i));
  }
  double *wlist[rank];
  SEXP ret = PROTECT(Rf_allocVector(VECSXP, rank));
  for (int i = 0; i < rank; i++) {
    SET_VECTOR_ELT(ret, i, Rf_allocVector(REALSXP, dims[i]));
    wlist[i] = REAL(VECTOR_ELT(ret, i));
  }

  SEXP Sdd = PROTECT(Rf_coerceVector(Sd, INTSXP));
  const int *dd = INTEGER_RO(Sdd);
  int threads = Rf_asInteger(Sthreads);
  for (int r = 0; r < rank; r++) {
    // The paper use index from 0 to n, i.e. n+1 knots in each dimension. Therefore dims[r]-1
    fhweights(dims[r] - 1, grid[r], dd[r], wlist[r], threads);
  }

  UNPROTECT(2);
  return ret;
}

double C_evalcheb(const double *cf, const double *x, const int *dims, const int rank) {
  if (rank == 0)
    return cf[0];
  // Otherwise, use the Clenshaw algorithm
  R_xlen_t siz = 1;
  const int newrank = rank - 1;
  const int N = dims[newrank];
  double x2 = 2.0 * x[newrank], bn1 = 0, bn2 = 0, bn = 0;

  // semantically unnecessary, just save a recursion level, it speeds up things.
  if (newrank == 0) {
    for (int i = N - 1; i > 0; i--) {
      bn2 = bn1;
      bn1 = bn;
      bn = x2 * bn1 - bn2 + cf[i];
    }
    return x[0] * bn - bn1 + cf[0];
  }
  for (int i = 0; i < newrank; i++)
    siz *= dims[i];
  for (int i = N - 1; i >= 0; i--) {
    bn2 = bn1;
    bn1 = bn;
    bn = x2 * bn1 - bn2 + C_evalcheb(&cf[i * siz], x, dims, newrank);
  }
  return bn - x[newrank] * bn1;
}

SEXP R_evalcheb(SEXP coef, SEXP inx, SEXP Rthreads, SEXP spare) {
  R_xlen_t siz = 1;
  const double *cf = REAL_RO(coef);
#ifdef _OPENMP
  int threads = Rf_asInteger(Rthreads);
#endif
  UNUSED(spare);
  SEXP dim = Rf_getAttrib(coef, R_DimSymbol);
  const int *dims = INTEGER_RO(dim);
  const int rank = LENGTH(dim);
  if (rank <= 0)
    Rf_error("rank must be positive");
  if (Rf_isMatrix(inx)) {
    if (rank != Rf_nrows(inx))
      Rf_error("coefficient rank(%d) must match number of rows(%d)", LENGTH(dim), Rf_nrows(inx));
  } else {
    if (rank != LENGTH(inx))
      Rf_error("coefficient rank(%d) does not match argument length(%d)", LENGTH(dim), LENGTH(inx));
  }
  // This shouldn't happen, but we check it anyway since we'll bomb if it's wrong
  for (int i = 0; i < rank; i++)
    siz *= dims[i];
  if (XLENGTH(coef) != siz)
    Rf_error("coefficient length(%lld) does not match data length(%lld)", (long long)XLENGTH(coef), (long long)siz);

  const double *xp = REAL_RO(inx);
  const int numvec = Rf_isMatrix(inx) ? Rf_ncols(inx) : 1;
  SEXP resvec = PROTECT(Rf_allocVector(REALSXP, numvec));
  double *out = REAL(resvec);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < numvec; i++) {
    out[i] = C_evalcheb(cf, xp + i * rank, dims, rank);
  }
  UNPROTECT(1);
  return resvec;
}

static double C_evalpolyh(const double *x, const double *knots, const double *weights, const double *lweights,
                          const int rank, const int nknots, const double k) {
  // Find squared distance to each knot
  int ki = (int)k;
  double wsum = 0.0;
  if (k < 0) {
    for (R_xlen_t i = 0; i < nknots; i++) {
      const double *kn = &knots[i * rank];
      double sqdist = 0.0;
      for (int j = 0; j < rank; j++)
        sqdist += (x[j] - kn[j]) * (x[j] - kn[j]);
      wsum += weights[i] * exp(k * sqdist);
    }
  } else if (ki % 2 == 1) {
    for (R_xlen_t i = 0; i < nknots; i++) {
      const double *kn = &knots[i * rank];
      double sqdist = 0.0;
      for (int j = 0; j < rank; j++)
        sqdist += (x[j] - kn[j]) * (x[j] - kn[j]);
      if (sqdist == 0)
        continue;
      wsum += weights[i] * R_pow_di(sqrt(sqdist), ki);
    }
  } else {
    for (R_xlen_t i = 0; i < nknots; i++) {
      const double *kn = &knots[i * rank];
      double sqdist = 0.0;
      for (int j = 0; j < rank; j++)
        sqdist += (x[j] - kn[j]) * (x[j] - kn[j]);
      if (sqdist == 0)
        continue;
      wsum += weights[i] * log(sqrt(sqdist)) * R_pow_di(sqrt(sqdist), ki);
    }
  }
  // then the linear part, constant is the first element
  wsum += lweights[0];
  for (int i = 0; i < rank; i++)
    wsum += lweights[i + 1] * x[i];
  return wsum;
}
SEXP R_evalpolyh(SEXP inx, SEXP Sknots, SEXP weights, SEXP lweights, SEXP Sk, SEXP Sthreads, SEXP spare) {
  const double *x = REAL_RO(inx);
  const double *knots = REAL_RO(Sknots);
  const double *w = REAL_RO(weights);
  const double *lw = REAL_RO(lweights);
  double k = Rf_asReal(Sk);
#ifdef _OPENMP
  int threads = Rf_asInteger(Sthreads);
#endif
  const int rank = Rf_nrows(Sknots);
  UNUSED(spare);
  if (LENGTH(lweights) != rank + 1)
    Rf_error("linear weights (%d) should be one longer than rank(%d)", LENGTH(lweights), rank);
  if (LENGTH(weights) != Rf_ncols(Sknots))
    Rf_error("number of weights (%d) should equal number of knots (%d)", LENGTH(weights), Rf_ncols(Sknots));
  if (Rf_isMatrix(inx)) {
    if (Rf_nrows(inx) != Rf_nrows(Sknots))
      Rf_error("Dimension of input (%d) must match dimension of interpolant (%d)", Rf_nrows(inx), Rf_nrows(Sknots));
  } else {
    if (Rf_length(inx) != Rf_nrows(Sknots))
      Rf_error("Dimension of input (%d) must match dimension of interpolant (%d)", Rf_length(inx), Rf_nrows(Sknots));
  }
  const int numvec = Rf_isMatrix(inx) ? Rf_ncols(inx) : 1;
  SEXP res = PROTECT(Rf_allocVector(REALSXP, numvec));
  double *out = REAL(res);
  const int nknots = Rf_ncols(Sknots);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < numvec; i++) {
    out[i] = C_evalpolyh(x + i * rank, knots, w, lw, rank, nknots, k);
  }
  UNPROTECT(1);
  return res;
}

// an R-free evalongrid. Recursive
void C_evalongrid(void (*fun)(double *x, double *y, int valuedim, void *ud), double *arg, const double **grid,
                  const int *dims, const int rank, const int valuedim, double *result, void *userdata) {
  int mrank = rank - 1;
  R_xlen_t stride = valuedim;

  if (rank == 0) {
    fun(arg, &result[0], valuedim, userdata);
    return;
  }

  for (int i = 0; i < mrank; i++)
    stride *= dims[i];
  for (int i = 0; i < dims[mrank]; i++) {
    arg[mrank] = grid[mrank][i];
    C_evalongrid(fun, arg, grid, dims, mrank, valuedim, &result[i * stride], userdata);
  }
}

static void C_call(double *x, double *y, const int valuedim, void *userdata) {
  // x is not needed: the argument vector is already part of the call in userdata
  UNUSED(x);
  PROTECT_INDEX ipx;
  SEXP fv = Rf_eval(*(SEXP *)userdata, R_BaseEnv);
  PROTECT_WITH_INDEX(fv, &ipx);
  if (!Rf_isNumeric(fv))
    Rf_error("fun must return a numeric value");
  if (LENGTH(fv) != valuedim)
    Rf_error("fun must return a value of length %d, not %d", valuedim, LENGTH(fv));
  REPROTECT(fv = Rf_coerceVector(fv, REALSXP), ipx);
  memcpy(y, REAL_RO(fv), valuedim * sizeof(double));
  UNPROTECT(1);
}

SEXP R_evalongrid(SEXP fun, SEXP sgrid) {
  if (!Rf_isFunction(fun))
    Rf_error("'fun' must be a function");
  const int rank = LENGTH(sgrid);
  const double *grid[rank];
  R_xlen_t len = 1;
  SEXP R_dim = PROTECT(Rf_allocVector(INTSXP, rank));
  int *dims = INTEGER(R_dim);
  for (int i = 0; i < rank; i++) {
    grid[i] = REAL_RO(VECTOR_ELT(sgrid, i));
    dims[i] = LENGTH(VECTOR_ELT(sgrid, i));
    if (dims[i] == 0)
      Rf_error("grid dimension %d has no points", i + 1);
    len *= dims[i];
  }

  // Create a call to fun
  SEXP R_fcall = PROTECT(Rf_lang2(fun, R_NilValue));
  // and an argument list to it
  SEXP R_arg = PROTECT(Rf_allocVector(REALSXP, rank));
  SETCADR(R_fcall, R_arg);

  // find the dimension of the result
  for (int i = 0; i < rank; i++)
    REAL(R_arg)[i] = grid[i][0];
  SEXP fv = Rf_eval(R_fcall, R_BaseEnv);
  if (!Rf_isNumeric(fv))
    Rf_error("fun must return a numeric value");
  const int valuedim = LENGTH(fv);

  // Now, it's just to fill in arg, and
  // do REAL(Rf_eval(R_fcall,NULL))[0]
  SEXP resvec = PROTECT(Rf_allocVector(REALSXP, len * valuedim));
  C_evalongrid(C_call, REAL(R_arg), grid, dims, rank, valuedim, REAL(resvec), (void *)&R_fcall);

  // Set the dim-attribute
  if (valuedim == 1) {
    Rf_setAttrib(resvec, R_DimSymbol, R_dim);
  } else {
    SEXP snewdim = PROTECT(Rf_allocVector(INTSXP, rank + 1));
    int *newdim = INTEGER(snewdim);
    newdim[0] = valuedim;
    memcpy(newdim + 1, dims, rank * sizeof(int));
    Rf_setAttrib(resvec, R_DimSymbol, snewdim);
    UNPROTECT(1);
  }
  UNPROTECT(4);
  return resvec;
}

// return index of the smallest element larger or equal to val
// binary search is slower than linear for small n, but we don't optimize that here
// small grids anyway run faster
// we never return 0, if val < arr[0] we return 1
// if val > arr[n-1], we return n-1
// this is suited for our particular application
static inline int binsearch(const int n, const double *arr, const double val) {
  int toosmall = 0, toolarge = n - 1;
  while (toosmall + 1 < toolarge) {
    const int test = (toosmall + toolarge) >> 1;
    if (arr[test] >= val) {
      toolarge = test;
    } else {
      toosmall = test;
    }
  }
  return toolarge;
}

double C_evalmlip(const int rank, const double *x, const double **grid, const int *dims, const double *values,
                  int blend) {

  double weight[rank];
  R_xlen_t valpos = 0;
  double ipval = 0;
  // Find first grid point which is larger or equal, compute
  // the weight, store it. As well as the linear position valpos in the grid
  // A binary search is faster if there are more grid points
  // but up to 10-20 linear search is usually equally fast or faster due to simplicity
  // Note that we never check whether x is within the grid. If it's outside we assume
  // the function continues linearly outwards.
  // Remaining optimizations:  We could have precomputed the stuff we divide by, if we should
  // do more than one vector in a row, this would save some time.  Likewise the stride (which is
  // integer multiplication, which is slow.)

  R_xlen_t stride = 1;
  for (int i = 0; i < rank; i++) {
    if (ISNAN(x[i]))
      return x[i];
    const double *gdvec = grid[i];
    int gp = binsearch(dims[i], gdvec, x[i]);
    weight[i] = (x[i] - gdvec[gp - 1]) / (gdvec[gp] - gdvec[gp - 1]);
    valpos += stride * gp;
    stride *= dims[i];
  }

  // loop over the corners of the box, sum values with weights
  double wsum = 0.0;
  for (int i = 0; i < (1 << rank); i++) {
    // i represents a corner. bit=1 if upper corner, 0 if lower corner.
    // We should find its weight
    // it's a product of the weights from each dimension
    R_xlen_t vpos = valpos;
    R_xlen_t stride = 1;
    double cw = 1;
    for (int g = 0; g < rank; g++) {
      if ((1 << g) & i) {
        cw *= blendfun(weight[g], blend);
      } else {
        cw *= blendfun(1 - weight[g], blend);
        vpos -= stride;
      }
      stride *= dims[g];
    }
    wsum += cw;
    ipval += cw * values[vpos];
  }
  return ipval / wsum;
}

// Then a multilinear approximation
SEXP R_evalmlip(SEXP sgrid, SEXP values, SEXP x, SEXP Rthreads, SEXP Sblend) {
  const int rank = LENGTH(sgrid);
  R_xlen_t gridsize = 1;
  int dims[rank];
#ifdef _OPENMP
  int threads = Rf_asInteger(Rthreads);
#endif
  const double *grid[rank];
  blend_type blend = BLEND_LINEAR;
  if (!Rf_isReal(values))
    Rf_error("values must be numeric");
  if (!Rf_isReal(x))
    Rf_error("argument x must be numeric");
  if (Rf_isMatrix(x) ? (Rf_nrows(x) != rank) : (LENGTH(x) != rank))
    Rf_error("grid has dimension %d, you supplied a length %d vector", rank, Rf_isMatrix(x) ? Rf_nrows(x) : LENGTH(x));

  // Make some pointers into the grid data
  for (int i = 0; i < rank; i++) {
    dims[i] = LENGTH(VECTOR_ELT(sgrid, i));
    grid[i] = REAL_RO(VECTOR_ELT(sgrid, i));
    gridsize *= dims[i];
  }

  if (XLENGTH(values) != gridsize)
    Rf_error("grid has size %lld, you supplied %lld values", (long long)gridsize, (long long)XLENGTH(values));
  if (!Rf_isNull(Sblend))
    blend = (blend_type)Rf_asInteger(Sblend);
  const int numvec = Rf_isMatrix(x) ? Rf_ncols(x) : 1;
  const double *xp = REAL_RO(x);
  SEXP resvec = PROTECT(Rf_allocVector(REALSXP, numvec));
  double *out = REAL(resvec);
  const double *val = REAL_RO(values);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < numvec; i++) {
    out[i] = C_evalmlip(rank, xp + i * rank, grid, dims, val, blend);
  }
  UNPROTECT(1);
  return resvec;
}

SEXP R_sqdiffs(SEXP x1, SEXP x2, SEXP Sthreads) {
  // each column in x1 should be subtracted from each column in x2,
  // the squared column sums should be returned.
  const int r1 = Rf_nrows(x1), c1 = Rf_ncols(x1), r2 = Rf_nrows(x2), c2 = Rf_ncols(x2);
  const R_xlen_t N = (R_xlen_t)c1 * c2;
#ifdef _OPENMP
  int threads = Rf_asInteger(Sthreads);
#endif
  SEXP res = PROTECT(Rf_allocVector(REALSXP, N));
  double *dres = REAL(res);
  double *np = dres;
  const double *x1v = REAL_RO(x1), *x2v = REAL_RO(x2);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (c1 > 1)
#endif
  for (R_xlen_t i = 0; i < c1; i++) {
    const double *x1p = x1v + i * r1;
    for (R_xlen_t j = 0; j < c2; j++) {
      const double *x2p = x2v + j * r2;
      double csum = 0.0;
      for (int k = 0; k < r1; k++) {
        csum += (x1p[k] - x2p[k]) * (x1p[k] - x2p[k]);
      }
      np[j + i * c2] = csum;
    }
  }
  SEXP snewdim = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(snewdim)[0] = c1;
  INTEGER(snewdim)[1] = c2;
  Rf_setAttrib(res, R_DimSymbol, snewdim);
  UNPROTECT(2);
  return res;
}

static double findsimplex(const double *x, const double *knots, const int *dtri, const double *lumats, const int *ipivs,
                          const double *bbox, bool epol, const double *val, const int dim, const int numsimplex,
                          int nknots, blend_type blend) {

  for (int i = 0; i < dim; i++)
    if (ISNAN(x[i]))
      return x[i];

  double vec[dim + 1];
  int N = dim + 1, one = 1, info;
  for (R_xlen_t simplex = 0; simplex < numsimplex; simplex++) {
    // x must be in the bounding box. Perhaps we should do binary search, or organize them hierarchically?
    const double *box = bbox + simplex * 2 * dim;
    bool bad = false;
    for (int i = 0; i < dim; i++) {
      if (x[i] < box[0] || x[i] > box[1]) {
        bad = true;
        break;
      }
      box += 2;
    }
    if (bad)
      continue;

    // bounding box matches. Transform to barycentric coordinates, check that they are positive.
    const double *lumat = lumats + simplex * N * N;
    const int *ipiv = ipivs + simplex * N;
    memcpy(vec, x, dim * sizeof(double));
    vec[dim] = 1.0;
    F77_CALL(dgetrs)("N", &N, &one, lumat, &N, ipiv, vec, &N, &info FCONE);
    for (int d = 0; d < N; d++)
      if (vec[d] < -1e-10) {
        bad = true;
        break;
      }
    if (bad)
      continue;

    // We found it
    const int *tri = dtri + simplex * N;
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
      const double w = blendfun(coord, blend);
      if (w == 0)
        continue;
      wsum += w;
      sum += w * val[tri[d] - 1];
    }
    return sum / wsum;
  }

  if (epol) {
    // Find the closest knot
    double mindist = INFINITY;
    int nearknot = 0;
    for (int k = 0; k < nknots; k++) {
      double dist = 0.0;
      for (int i = 0; i < dim; i++)
        dist += (x[i] - knots[(R_xlen_t)k * dim + i]) * (x[i] - knots[(R_xlen_t)k * dim + i]);
      if (dist < mindist) {
        mindist = dist;
        nearknot = k;
      }
    }

    // Now, find a simplex with this knot
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
    const int *tri = dtri + nearest * (dim + 1);
    memcpy(vec, x, dim * sizeof(double));
    vec[dim] = 1.0;
    F77_CALL(dgetrs)("N", &N, &one, lumat, &N, ipiv, vec, &N, &info FCONE);
    double sum = 0;
    for (int d = 0; d <= dim; d++)
      sum += vec[d] * val[tri[d] - 1];
    return sum;
  }
  return NA_REAL;
}

SEXP R_evalsl(SEXP Sx, SEXP Sknots, SEXP Sdtri, SEXP adata, SEXP Sval, SEXP extrapolate, SEXP Sthreads, SEXP Sblend) {
  // Loop over the triangulation
  const int *dtri = INTEGER_RO(Sdtri);
  const int numsimplex = Rf_ncols(Sdtri);
  const double *lumats = REAL_RO(VECTOR_ELT(adata, 0));
  const int *pivots = INTEGER_RO(VECTOR_ELT(adata, 1));
  const double *bboxmat = REAL_RO(VECTOR_ELT(adata, 2));
  const double *x = REAL_RO(Sx);
  const int dim = Rf_nrows(Sknots);
  const int numvec = Rf_ncols(Sx);
  if (Rf_nrows(Sx) != dim)
    Rf_error("Dimension of input (%d) must match dimension of interpolant (%d)", Rf_nrows(Sx), dim);
  const double *knots = REAL_RO(Sknots);
  const double *val = REAL_RO(Sval);
#ifdef _OPENMP
  int threads = Rf_asInteger(Sthreads);
#endif
  SEXP ret = PROTECT(Rf_allocVector(REALSXP, Rf_ncols(Sx)));
  double *resvec = REAL(ret);
  bool epol = Rf_asLogical(extrapolate);
  blend_type blend = BLEND_LINEAR;
  if (!Rf_isNull(Sblend))
    blend = (blend_type)Rf_asInteger(Sblend);
  const int nknots = Rf_ncols(Sknots);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(guided) if (threads > 1 && numvec > 1)
#endif
  for (R_xlen_t i = 0; i < numvec; i++) {
    resvec[i] =
        findsimplex(x + i * dim, knots, dtri, lumats, pivots, bboxmat, epol, val, dim, numsimplex, nknots, blend);
  }
  UNPROTECT(1);
  return ret;
}

SEXP R_analyzesimplex(SEXP Sdtri, SEXP Sknots, SEXP Sthreads) {
#ifdef _OPENMP
  int threads = Rf_asInteger(Sthreads);
#endif
  const int *dtri = INTEGER_RO(Sdtri);
  const double *knots = REAL_RO(Sknots);
  const int dim = Rf_nrows(Sknots);
  const int numsimplex = Rf_ncols(Sdtri);
  int N = dim + 1;
  SEXP retlist = PROTECT(Rf_allocVector(VECSXP, 3));
  // Allocate everything before parallel for
  SET_VECTOR_ELT(retlist, 0, Rf_allocVector(REALSXP, (R_xlen_t)numsimplex * N * N));
  SET_VECTOR_ELT(retlist, 1, Rf_allocVector(INTSXP, (R_xlen_t)numsimplex * N));
  SET_VECTOR_ELT(retlist, 2, Rf_allocMatrix(REALSXP, 2 * dim, numsimplex));
  double *lumats = REAL(VECTOR_ELT(retlist, 0));
  int *pivots = INTEGER(VECTOR_ELT(retlist, 1));
  double *bboxmat = REAL(VECTOR_ELT(retlist, 2));

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1 && numsimplex > 1)
#endif
  for (R_xlen_t simplex = 0; simplex < numsimplex; simplex++) {
    const int *tri = dtri + simplex * (dim + 1);
    int *ipiv = pivots + simplex * N;
    for (int j = 0; j < N; j++)
      ipiv[j] = j + 1;
    double *lumat = lumats + simplex * N * N;
    int info;
    double *m = bboxmat + simplex * 2 * dim;
    for (int j = 0; j < dim; j++) {
      m[2 * j] = INFINITY;
      m[2 * j + 1] = -INFINITY;
    }
    // Copy all vertices into mat, with a final row of ones to LU-factorize
    // Preparation for transforming into barycentric coordinates
    for (int v = 0; v <= dim; v++) {
      const double *knot = knots + (tri[v] - 1) * dim;
      for (int j = 0; j < dim; j++) {
        lumat[v * (dim + 1) + j] = knot[j];
        if (knot[j] < m[2 * j])
          m[2 * j] = knot[j];
        if (knot[j] > m[2 * j + 1])
          m[2 * j + 1] = knot[j];
      }
      lumat[v * (dim + 1) + dim] = 1.0;
    }
    F77_CALL(dgetrf)(&N, &N, lumat, &N, ipiv, &info);
  }

  UNPROTECT(1);
  return retlist;
}

SEXP R_havefftw(void) {
#ifdef HAVE_FFTW
  return Rf_ScalarLogical(TRUE);
#else
  return Rf_ScalarLogical(FALSE);
#endif
}

// inplace to save memory, do repeated multiplication instead of pow(). It's faster.
SEXP R_phifunc(SEXP Sx, SEXP Sk, SEXP Sthreads) {
  double k = Rf_asReal(Sk);
  double *x = REAL(Sx);
  const R_xlen_t n = XLENGTH(Sx);
#ifdef _OPENMP
  int threads = Rf_asInteger(Sthreads);
#endif
  if (n == 1 && x[0] == 0.0)
    return Rf_ScalarReal((k < 0) ? 1 : 0);

  double *y;
  SEXP res;
  if (MAYBE_REFERENCED(Sx)) {
    res = PROTECT(Rf_allocVector(REALSXP, n));
    y = REAL(res);
  } else {
    res = Sx;
    y = x;
  }

  if (k < 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
#endif
    for (R_xlen_t i = 0; i < n; i++)
      y[i] = exp(k * x[i]);
  } else {
    int ki = Rf_asInteger(Sk);
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
  return res;
}
#ifndef HAVE_ALGLIB
SEXP R_makerbf(SEXP a, SEXP b, SEXP c, SEXP d) {
  UNUSED(a);
  UNUSED(b);
  UNUSED(c);
  UNUSED(d);
  Rf_error("ALGLIB not available");
  return R_NilValue;
}
SEXP R_evalrbf(SEXP a, SEXP b, SEXP c) {
  UNUSED(a);
  UNUSED(b);
  UNUSED(c);
  Rf_error("ALGLIB not available");
  return R_NilValue;
}
SEXP R_havealglib(void) { return Rf_ScalarLogical(FALSE); }

#endif
