#include "chebpol.h"

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

static const R_CallMethodDef CallEntries[] = {
    {"c_eval_cheb", (DL_FUNC)&c_eval_cheb, 4},
    {"c_cheb_coef", (DL_FUNC)&c_cheb_coef, 3},
    {"c_fh", (DL_FUNC)&c_fh, 6},
    {"c_fh_weights", (DL_FUNC)&c_fh_weights, 3},
    {"c_eval_mlip", (DL_FUNC)&c_eval_mlip, 5},
    {"c_eval_on_grid", (DL_FUNC)&c_eval_on_grid, 2},
    {"c_have_fftw", (DL_FUNC)&c_have_fftw, 0},
    {"c_sq_diffs", (DL_FUNC)&c_sq_diffs, 3},
    {"c_phi_func", (DL_FUNC)&c_phi_func, 3},
    {"c_make_rbf", (DL_FUNC)&c_make_rbf, 4},
    {"c_eval_rbf", (DL_FUNC)&c_eval_rbf, 3},
    {"c_have_alglib", (DL_FUNC)&c_have_alglib, 0},
    {"c_eval_polyh", (DL_FUNC)&c_eval_polyh, 7},
    {"c_eval_sl", (DL_FUNC)&c_eval_sl, 8},
    {"c_eval_stalker", (DL_FUNC)&c_eval_stalker, 5},
    {"c_make_hyp", (DL_FUNC)&c_make_hyp, 2},
    {"c_eval_hyp", (DL_FUNC)&c_eval_hyp, 4},
    {"c_make_stalk", (DL_FUNC)&c_make_stalk, 2},
    {"c_eval_stalk", (DL_FUNC)&c_eval_stalk, 4},
    {"c_analyze_simplex", (DL_FUNC)&c_analyze_simplex, 3},
    {"c_have_gsl", (DL_FUNC)&c_have_gsl, 0},
    {NULL, NULL, 0},
};

void attribute_visible R_init_chebpol(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
  R_RegisterCCallable("chebpol", "evalcheb", (DL_FUNC)&eval_cheb);
  R_RegisterCCallable("chebpol", "evalmlip", (DL_FUNC)&eval_mlip);
  R_RegisterCCallable("chebpol", "evalongrid", (DL_FUNC)&eval_on_grid);
}

#ifdef HAVE_FFTW
void attribute_visible R_unload_chebpol(DllInfo *dll) {
  (void)dll;
  fftw_cleanup();
}
#endif
