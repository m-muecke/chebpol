#include "chebpol.h"

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

static const R_CallMethodDef CallEntries[] = {
    {"c_eval_chebyshev", (DL_FUNC)&c_eval_chebyshev, 3},
    {"c_chebyshev_coef", (DL_FUNC)&c_chebyshev_coef, 3},
    {"c_eval_fh", (DL_FUNC)&c_eval_fh, 5},
    {"c_fh_weights", (DL_FUNC)&c_fh_weights, 3},
    {"c_eval_multilinear", (DL_FUNC)&c_eval_multilinear, 5},
    {"c_eval_on_grid", (DL_FUNC)&c_eval_on_grid, 2},
    {"c_has_fftw", (DL_FUNC)&c_has_fftw, 0},
    {"c_sq_dists", (DL_FUNC)&c_sq_dists, 3},
    {"c_rbf_kernel", (DL_FUNC)&c_rbf_kernel, 3},
    {"c_eval_polyharmonic", (DL_FUNC)&c_eval_polyharmonic, 6},
    {"c_eval_simplexlinear", (DL_FUNC)&c_eval_simplexlinear, 8},
    {"c_eval_oldstalker", (DL_FUNC)&c_eval_oldstalker, 5},
    {"c_make_hstalker", (DL_FUNC)&c_make_hstalker, 2},
    {"c_eval_hstalker", (DL_FUNC)&c_eval_hstalker, 4},
    {"c_make_stalker", (DL_FUNC)&c_make_stalker, 2},
    {"c_eval_stalker", (DL_FUNC)&c_eval_stalker, 4},
    {"c_make_simplexlinear", (DL_FUNC)&c_make_simplexlinear, 3},
    {"c_has_gsl", (DL_FUNC)&c_has_gsl, 0},
    {NULL, NULL, 0},
};

void attribute_visible R_init_chebpol(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

#ifdef HAVE_FFTW
void attribute_visible R_unload_chebpol(DllInfo *dll) {
  (void)dll;
  fftw_cleanup();
}
#endif
