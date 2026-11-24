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
// On Windows configure is not run; Makevars.win sets the HAVE_ macros instead.
#ifndef _WIN32
#include "config.h"
#endif

#ifndef FCONE
#define FCONE
#endif

/* .Call entry points, registered in init.c */
SEXP R_evalcheb(SEXP, SEXP, SEXP, SEXP);
SEXP R_chebcoef(SEXP, SEXP, SEXP);
SEXP R_FH(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_fhweights(SEXP, SEXP, SEXP);
SEXP R_evalmlip(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_mlippred(SEXP, SEXP);
SEXP R_evalongrid(SEXP, SEXP);
SEXP R_havefftw(void);
SEXP R_sqdiffs(SEXP, SEXP, SEXP);
SEXP R_phifunc(SEXP, SEXP, SEXP);
SEXP R_makerbf(SEXP, SEXP, SEXP, SEXP);
SEXP R_evalrbf(SEXP, SEXP, SEXP);
SEXP R_havealglib(void);
SEXP R_evalpolyh(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_evalsl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_evalstalker(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_makehyp(SEXP, SEXP);
SEXP R_evalhyp(SEXP, SEXP, SEXP, SEXP);
SEXP R_makestalk(SEXP, SEXP);
SEXP R_evalstalk(SEXP, SEXP, SEXP, SEXP);
SEXP R_analyzesimplex(SEXP, SEXP, SEXP);
SEXP havegsl(void);

/* C entry points made available to other packages in init.c */
double C_evalcheb(const double *cf, const double *x, const int *dims, const int rank);
double C_evalmlip(const int rank, const double *x, const double **grid, const int *dims, const double *values,
                  int blend);
void C_evalongrid(void (*fun)(double *x, double *y, int valuedim, void *ud), double *arg, const double **grid,
                  const int *dims, const int rank, const int valuedim, double *result, void *userdata);

static R_INLINE double blendfun(double w, int blend) {
  switch (blend) {
  case 0:
    // Linear, nothing to do here
    break;
  case 1:
    w = w < 0.5 ? 0.5 * exp(2 - 1 / w) : 1 - 0.5 * exp(2 - 1 / (1 - w)); // smooth sigmoid blending
    break;
  case 2:
    w = w < 0.5 ? 0.5 * exp(4 - 1 / (w * w)) : 1 - 0.5 * exp(4 - 1 / ((1 - w) * (1 - w))); // parodic sigmoid blending
    break;
  case 3:
    w = (-2 * w + 3) * w * w; // cubic blending
    break;
  case 4:
    w = (w < 0.5) ? 0 : 1; // discontinuous blending
    break;
  case 5:
    w = 0.5; // mean blending
  }
  return w;
}
#endif
