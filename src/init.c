#include "chebpol.h"

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

static const R_CallMethodDef CallEntries[] = {
    {"evalcheb", (DL_FUNC)&R_evalcheb, 4},
    {"chebcoef", (DL_FUNC)&R_chebcoef, 3},
    {"FH", (DL_FUNC)&R_FH, 6},
    {"FHweights", (DL_FUNC)&R_fhweights, 3},
    {"evalmlip", (DL_FUNC)&R_evalmlip, 5},
    {"evalongrid", (DL_FUNC)&R_evalongrid, 2},
    {"havefftw", (DL_FUNC)&R_havefftw, 0},
    {"sqdiffs", (DL_FUNC)&R_sqdiffs, 3},
    {"phifunc", (DL_FUNC)&R_phifunc, 3},
    {"makerbf", (DL_FUNC)&R_makerbf, 4},
    {"evalrbf", (DL_FUNC)&R_evalrbf, 3},
    {"havealglib", (DL_FUNC)&R_havealglib, 0},
    {"evalpolyh", (DL_FUNC)&R_evalpolyh, 7},
    {"evalsl", (DL_FUNC)&R_evalsl, 8},
    {"evalstalker", (DL_FUNC)&R_evalstalker, 5},
    {"makehyp", (DL_FUNC)&R_makehyp, 2},
    {"evalhyp", (DL_FUNC)&R_evalhyp, 4},
    {"makestalk", (DL_FUNC)&R_makestalk, 2},
    {"evalstalk", (DL_FUNC)&R_evalstalk, 4},
    {"analyzesimplex", (DL_FUNC)&R_analyzesimplex, 3},
    {"havegsl", (DL_FUNC)&havegsl, 0},
    {NULL, NULL, 0},
};

void attribute_visible R_init_chebpol(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_RegisterCCallable("chebpol", "evalcheb", (DL_FUNC)&C_evalcheb);
  R_RegisterCCallable("chebpol", "evalmlip", (DL_FUNC)&C_evalmlip);
  R_RegisterCCallable("chebpol", "evalongrid", (DL_FUNC)&C_evalongrid);
}

#ifdef HAVE_FFTW
void attribute_visible R_unload_chebpol(DllInfo *dll) {
  (void)dll;
  fftw_cleanup();
}
#endif
