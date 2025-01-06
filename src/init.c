#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(favas)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mace )(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(set_alpha)(void *);
extern void F77_NAME(set_big  )(void *);
extern void F77_NAME(set_span )(void *);
extern void F77_NAME(set_sml  )(void *);
extern void F77_NAME(set_eps  )(void *);
extern void F77_NAME(set_spans)(void *);
extern void F77_NAME(set_maxit)(void *); 
extern void F77_NAME(set_nterm)(void *); 


static const R_FortranMethodDef FortranEntries[] =
{
  {"favas",     (DL_FUNC) &F77_NAME(favas),     16},
  {"mace",      (DL_FUNC) &F77_NAME(mace),      14},
  {"set_alpha", (DL_FUNC) &F77_NAME(set_alpha),  1},
  {"set_big",   (DL_FUNC) &F77_NAME(set_big),    1},
  {"set_span",  (DL_FUNC) &F77_NAME(set_span),   1},
  {"set_sml",   (DL_FUNC) &F77_NAME(set_sml),    1},
  {"set_eps",   (DL_FUNC) &F77_NAME(set_eps),    1},
  {"set_spans", (DL_FUNC) &F77_NAME(set_spans),  1},
  {"set_maxit", (DL_FUNC) &F77_NAME(set_maxit),  1}, 
  {"set_nterm", (DL_FUNC) &F77_NAME(set_nterm),  1}, 
  {NULL, NULL, 0}
};

void R_init_acepack(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
