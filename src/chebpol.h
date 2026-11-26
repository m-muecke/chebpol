#ifndef CHEBPOL_H
#define CHEBPOL_H
#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <float.h>
#include <stdbool.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <R_ext/Visibility.h>

#ifndef FCONE
#define FCONE
#endif

// .Call entry points, registered in init.c
SEXP c_eval_cheb(SEXP, SEXP, SEXP, SEXP);
SEXP c_cheb_coef(SEXP, SEXP, SEXP);
SEXP c_fh(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP c_fh_weights(SEXP, SEXP, SEXP);
SEXP c_eval_mlip(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP c_eval_on_grid(SEXP, SEXP);
SEXP c_has_fftw(void);
SEXP c_sq_diffs(SEXP, SEXP, SEXP);
SEXP c_phi_func(SEXP, SEXP, SEXP);
SEXP c_make_rbf(SEXP, SEXP, SEXP, SEXP);
SEXP c_eval_rbf(SEXP, SEXP, SEXP);
SEXP c_has_alglib(void);
SEXP c_eval_polyh(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP c_eval_sl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP c_eval_oldstalker(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP c_make_hstalker(SEXP, SEXP);
SEXP c_eval_hstalker(SEXP, SEXP, SEXP, SEXP);
SEXP c_make_stalker(SEXP, SEXP);
SEXP c_eval_stalker(SEXP, SEXP, SEXP, SEXP);
SEXP c_analyze_simplex(SEXP, SEXP, SEXP);
SEXP c_has_gsl(void);

// C entry points made available to other packages in init.c
double eval_cheb(const double *cf, const double *x, const int *dims, const int rank);
double eval_mlip(const int rank, const double *x, const double **grid, const int *dims, const double *values,
                 int blend);
void eval_on_grid(void (*fun)(double *x, double *y, int valuedim, void *ud), double *arg, const double **grid,
                  const int *dims, const int rank, const int valuedim, double *result, void *userdata);

typedef enum {
  BLEND_LINEAR = 0,
  BLEND_SIGMOID = 1,
  BLEND_PARODIC = 2,
  BLEND_CUBIC = 3,
  BLEND_SQUARE = 4,
  BLEND_MEAN = 5
} blend_type;

// Index i of the grid interval [xt[i], xt[i + 1]) that holds x, clamped to 0 .. n - 2.
static inline int find_interval(const double *xt, const int n, const double x) {
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (x >= xt[mid])
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

static inline double blend_fun(double w, blend_type blend) {
  // only the linear blend is defined outside [0, 1], where it extrapolates
  if (blend != BLEND_LINEAR)
    w = w < 0 ? 0 : (w > 1 ? 1 : w);
  switch (blend) {
  case BLEND_LINEAR:
    break;
  case BLEND_SIGMOID:
    w = w < 0.5 ? 0.5 * exp(2 - 1 / w) : 1 - 0.5 * exp(2 - 1 / (1 - w)); // smooth sigmoid blending
    break;
  case BLEND_PARODIC:
    w = w < 0.5 ? 0.5 * exp(4 - 1 / (w * w)) : 1 - 0.5 * exp(4 - 1 / ((1 - w) * (1 - w))); // parodic sigmoid blending
    break;
  case BLEND_CUBIC:
    w = (-2 * w + 3) * w * w; // cubic blending
    break;
  case BLEND_SQUARE:
    w = (w < 0.5) ? 0 : 1; // discontinuous blending
    break;
  case BLEND_MEAN:
    w = 0.5; // mean blending
  }
  return w;
}
#endif
