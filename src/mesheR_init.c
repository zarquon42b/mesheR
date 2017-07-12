#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP angcheck(SEXP, SEXP, SEXP);
extern SEXP areafun(SEXP);
extern SEXP displaceGauss(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP multisolve3(SEXP, SEXP);
extern SEXP smoothField(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP trianvol(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"angcheck",      (DL_FUNC) &angcheck,       3},
    {"areafun",       (DL_FUNC) &areafun,        1},
    {"displaceGauss", (DL_FUNC) &displaceGauss, 17},
    {"multisolve3",   (DL_FUNC) &multisolve3,    2},
    {"smoothField",   (DL_FUNC) &smoothField,   10},
    {"trianvol",      (DL_FUNC) &trianvol,       3},
    {NULL, NULL, 0}
};

void R_init_mesheR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
